# export_triad_fields.py
# Run with:
#   abaqus python export_triad_fields.py --odb m1.odb --out triad_fields.inp
#
# Optional:
#   --step Step-1
#   --nset ALL_NODES
#   --axis-point 0 0 0
#   --axis-dir 0 0 1
#   --m-field NT11
#
# Output fields:
#   VARIABLE=1  -> m
#   VARIABLE=2  -> e_r_x
#   VARIABLE=3  -> e_r_y
#   VARIABLE=4  -> e_r_z
#   VARIABLE=5  -> e_theta_x
#   VARIABLE=6  -> e_theta_y
#   VARIABLE=7  -> e_theta_z
#   VARIABLE=8  -> e_z_x
#   VARIABLE=9  -> e_z_y
#   VARIABLE=10 -> e_z_z

from __future__ import print_function

import argparse
import math
from collections import OrderedDict

from odbAccess import openOdb
from abaqusConstants import NODAL


# ----------------------------------------------------------------------
# Small vector helpers used to construct the geometric triad.
# ----------------------------------------------------------------------

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def norm(a):
    return math.sqrt(max(dot(a, a), 0.0))


def sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])


def scale(scalar, vector):
    return (scalar*vector[0], scalar*vector[1], scalar*vector[2])


def cross(a, b):
    return (
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    )


def unit(vector, fallback=None):
    vector_norm = norm(vector)
    if vector_norm < 1.0e-14:
        if fallback is None:
            raise ValueError("Cannot normalize near-zero vector.")
        return fallback
    return (
        vector[0]/vector_norm,
        vector[1]/vector_norm,
        vector[2]/vector_norm
    )



# ----------------------------------------------------------------------
# Command-line interface.
# ----------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--odb", required=True, help="Path to Model 1 ODB")
    parser.add_argument("--out", default="triad_fields.inp",
                        help="Output include file")
    parser.add_argument("--step", default=None,
                        help="Step name; default is last step")
    parser.add_argument("--nset", default="ALL_NODES",
                        help="Node set to export")
    parser.add_argument(
        "--axis-point", nargs=3, type=float, default=[0.0, 0.0, 0.0],
        metavar=("X0", "Y0", "Z0"),
        help="A point on the tube axis"
    )
    parser.add_argument(
        "--axis-dir", nargs=3, type=float, default=[0.0, 0.0, 1.0],
        metavar=("DX", "DY", "DZ"),
        help="Tube axis direction"
    )
    parser.add_argument(
        "--m-field", default="NT11",
        help="Field output name to use for fibrosis scalar m"
    )
    return parser.parse_args()



# ----------------------------------------------------------------------
# ODB access helpers.
# ----------------------------------------------------------------------

def get_last_step(odb):
    step_names = list(odb.steps.keys())
    if not step_names:
        raise RuntimeError("ODB has no steps.")
    return odb.steps[step_names[-1]]


def get_nset_instance_and_labels(root_assembly, nset_name):
    if nset_name in root_assembly.nodeSets:
        node_set = root_assembly.nodeSets[nset_name]
        labels_by_instance = OrderedDict()
        for node in node_set.nodes:
            instance_name = node.instanceName
            labels_by_instance.setdefault(instance_name, set()).add(
                node.label
            )
        return labels_by_instance

    labels_by_instance = OrderedDict()
    found = False
    for instance_name, instance in root_assembly.instances.items():
        if nset_name in instance.nodeSets:
            found = True
            node_set = instance.nodeSets[nset_name]
            labels_by_instance.setdefault(instance_name, set())
            for node in node_set.nodes:
                labels_by_instance[instance_name].add(node.label)

    if not found:
        raise RuntimeError("Node set '{}' not found.".format(nset_name))

    return labels_by_instance


def collect_deformed_coords(root_assembly, frame, labels_by_instance):
    deformed_coords = {}
    displacement_field = frame.fieldOutputs["U"]

    for instance_name, labels in labels_by_instance.items():
        instance = root_assembly.instances[instance_name]
        subset = displacement_field.getSubset(region=instance,
                                              position=NODAL)

        displacement_by_label = {}
        for value in subset.values:
            if value.nodeLabel in labels:
                displacement_by_label[value.nodeLabel] = tuple(value.data)

        for node in instance.nodes:
            if node.label not in labels:
                continue
            reference_position = tuple(node.coordinates)
            displacement = displacement_by_label.get(
                node.label, (0.0, 0.0, 0.0)
            )
            current_position = add(reference_position, displacement)
            deformed_coords[node.label] = current_position

    return deformed_coords


def collect_scalar_field(root_assembly, frame, labels_by_instance,
                         field_name):
    if field_name not in frame.fieldOutputs:
        raise RuntimeError("Field '{}' not found in last frame.".format(
            field_name
        ))

    scalar_by_label = {}
    field_output = frame.fieldOutputs[field_name]

    for instance_name, labels in labels_by_instance.items():
        instance = root_assembly.instances[instance_name]
        subset = field_output.getSubset(region=instance, position=NODAL)

        for value in subset.values:
            if value.nodeLabel in labels:
                data = value.data
                if isinstance(data, (tuple, list)):
                    if len(data) != 1:
                        raise RuntimeError(
                            "Field '{}' is not scalar at node {}.".format(
                                field_name, value.nodeLabel
                            )
                        )
                    data = data[0]
                scalar_by_label[value.nodeLabel] = float(data)

    return scalar_by_label



# ----------------------------------------------------------------------
# Build e_r, e_theta, e_z from deformed coordinates.
# ----------------------------------------------------------------------

def geometric_triad(current_position, axis_point, e_z):
    # Project point onto axis, then define radial direction from axis to point.
    position_from_axis_point = sub(current_position, axis_point)
    axial_component = scale(dot(position_from_axis_point, e_z), e_z)
    radial_vector = sub(position_from_axis_point, axial_component)

    # If the node sits on the tube axis, radial_vector collapses to ~0 and
    # the radial direction is undefined. To keep the triad well-defined, pick
    # an arbitrary reference vector and strip off its axial part, leaving a
    # valid radial direction in the plane perpendicular to e_z.
    if norm(radial_vector) < 1.0e-12:
        trial_vector = (1.0, 0.0, 0.0)
        # Use x as that reference, unless x is nearly parallel to the axis,
        # in which case stripping the axial part would leave almost nothing
        # (ill-conditioned). |dot| > 0.9 means x is within ~26 deg of e_z, so
        # switch to y, which is then safely perpendicular to the axis.
        if abs(dot(trial_vector, e_z)) > 0.9:
            trial_vector = (0.0, 1.0, 0.0)
        radial_vector = sub(
            trial_vector,
            scale(dot(trial_vector, e_z), e_z)
        )

    e_r = unit(radial_vector)
    e_theta = unit(cross(e_z, e_r))
    e_r = unit(cross(e_theta, e_z))
    return e_r, e_theta, e_z



# ----------------------------------------------------------------------
# Abaqus predefined-field writer.
# ----------------------------------------------------------------------

def write_field_block(file_handle, variable_id, data_by_label):
    file_handle.write("*INITIAL CONDITIONS, TYPE=FIELD, VARIABLE={}\n"
                      .format(variable_id))
    for label in sorted(data_by_label):
        file_handle.write("{}, {:.16g}\n".format(label,
                                                  data_by_label[label]))



# ----------------------------------------------------------------------
# Main export workflow.
# ----------------------------------------------------------------------

def main():
    args = parse_args()

    odb = openOdb(args.odb, readOnly=True)
    try:
        step = odb.steps[args.step] if args.step else get_last_step(odb)
        if len(step.frames) == 0:
            raise RuntimeError("Selected step has no frames.")
        frame = step.frames[-1]

        axis_point = tuple(args.axis_point)
        e_z_axis = unit(tuple(args.axis_dir), fallback=(0.0, 0.0, 1.0))

        labels_by_instance = get_nset_instance_and_labels(
            odb.rootAssembly, args.nset
        )
        deformed_coords = collect_deformed_coords(
            odb.rootAssembly, frame, labels_by_instance
        )
        m_values = collect_scalar_field(
            odb.rootAssembly, frame, labels_by_instance, args.m_field
        )

        missing_m_labels = [
            label for label in deformed_coords if label not in m_values
        ]
        if missing_m_labels:
            raise RuntimeError(
                "Missing '{}' values for {} nodes. Example label: {}"
                .format(args.m_field, len(missing_m_labels),
                        missing_m_labels[0])
            )

        e_r_x = {}
        e_r_y = {}
        e_r_z = {}
        e_theta_x = {}
        e_theta_y = {}
        e_theta_z = {}
        e_z_x = {}
        e_z_y = {}
        e_z_z = {}

        for label, current_position in deformed_coords.items():
            e_r, e_theta, e_z = geometric_triad(current_position,
                                                axis_point, e_z_axis)
            e_r_x[label], e_r_y[label], e_r_z[label] = e_r
            e_theta_x[label], e_theta_y[label], e_theta_z[label] = e_theta
            e_z_x[label], e_z_y[label], e_z_z[label] = e_z

        with open(args.out, "w") as output_file:
            write_field_block(output_file, 1, m_values)
            write_field_block(output_file, 2, e_r_x)
            write_field_block(output_file, 3, e_r_y)
            write_field_block(output_file, 4, e_r_z)
            write_field_block(output_file, 5, e_theta_x)
            write_field_block(output_file, 6, e_theta_y)
            write_field_block(output_file, 7, e_theta_z)
            write_field_block(output_file, 8, e_z_x)
            write_field_block(output_file, 9, e_z_y)
            write_field_block(output_file, 10, e_z_z)

        print("Wrote {}".format(args.out))
        print("Fields:")
        print("  1 = m")
        print("  2:4 = e_r")
        print("  5:7 = e_theta")
        print("  8:10 = e_z")

    finally:
        odb.close()


if __name__ == "__main__":
    main()
