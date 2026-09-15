# motility_from_cross_sections.py
# Run:
#   abaqus viewer noGUI=motility_from_cross_sections.py
#
# Slice-integration method:
# - read reference node coordinates and INNER_SURFACE from mesh.inp
# - use reference z levels in INNER_SURFACE as cross-section buckets
# - read COORD field output frame-by-frame from the ODB
# - compute cross-sectional area in each xy slice using the cylinder axis
# - integrate areas along z with the trapezoidal rule
#
# Output:
#   motility_volume_cross_sections.csv

from odbAccess import openOdb
import csv
import math

INP_PATH = 'mesh.inp'
ODB_PATH = 'm2.odb'
STEP_NAME = 'Step-1'
NODESET_NAME = 'INNER_SURFACE'
CSV_OUT = 'motility_volume_cross_sections.csv'

Z_TOL = 1.0e-8

# Nodes whose reference z coordinates differ by less than this tolerance
# are treated as belonging to the same cross-section.


# ----------------------------------------------------------------------
# Basic helpers
# ----------------------------------------------------------------------

def mean(values):
    value_sum = 0.0
    index = 0
    while index < len(values):
        value_sum += values[index]
        index += 1
    return value_sum / float(len(values))


# ----------------------------------------------------------------------
# Parse mesh.inp
# ----------------------------------------------------------------------

def parse_inp_nodes_and_nset(input_path, node_set_name):
    reference_nodes = {}
    labels_in_node_set = []

    with open(input_path, 'r') as input_file:
        lines = input_file.readlines()

    line_index = 0
    while line_index < len(lines):
        raw_line = lines[line_index].strip()
        lower_line = raw_line.lower()

        is_node_block = (
            lower_line.startswith('*node')
            and not lower_line.startswith('*node output')
            and not lower_line.startswith('*node print')
        )
        if is_node_block:
            line_index += 1
            while (line_index < len(lines)
                   and not lines[line_index].strip().startswith('*')):
                fields = [
                    field.strip()
                    for field in lines[line_index].split(',')
                    if field.strip()
                ]
                if len(fields) >= 4:
                    node_label = int(fields[0])
                    reference_nodes[node_label] = (
                        float(fields[1]),
                        float(fields[2]),
                        float(fields[3])
                    )
                line_index += 1
            continue

        is_requested_nset = (
            lower_line.startswith('*nset')
            and ('nset=' + node_set_name.lower()) in lower_line
        )
        if is_requested_nset:
            generate = ('generate' in lower_line)
            line_index += 1
            while (line_index < len(lines)
                   and not lines[line_index].strip().startswith('*')):
                fields = [
                    field.strip()
                    for field in lines[line_index].split(',')
                    if field.strip()
                ]
                if generate and len(fields) >= 3:
                    start_label = int(fields[0])
                    end_label = int(fields[1])
                    label_step = int(fields[2])
                    node_label = start_label
                    while node_label <= end_label:
                        labels_in_node_set.append(node_label)
                        node_label += label_step
                else:
                    field_index = 0
                    while field_index < len(fields):
                        labels_in_node_set.append(int(fields[field_index]))
                        field_index += 1
                line_index += 1
            continue

        line_index += 1

    labels_in_node_set = sorted(set(labels_in_node_set))

    if len(reference_nodes) == 0:
        raise RuntimeError('No nodes parsed from mesh.inp')
    if len(labels_in_node_set) == 0:
        raise RuntimeError('Node set %s not found in mesh.inp'
                           % node_set_name)

    return reference_nodes, labels_in_node_set


# ----------------------------------------------------------------------
# Build cross-section buckets from reference z
# ----------------------------------------------------------------------


# The z buckets define the axial integration stations for the lumen
# volume calculation.  Node ordering within each bucket is angular.
def build_z_buckets(reference_nodes, inner_surface_labels):
    groups_by_z_key = {}

    index = 0
    while index < len(inner_surface_labels):
        node_label = inner_surface_labels[index]
        z_value = reference_nodes[node_label][2]
        z_key = int(round(z_value / Z_TOL))
        if z_key not in groups_by_z_key:
            groups_by_z_key[z_key] = []
        groups_by_z_key[z_key].append(node_label)
        index += 1

    z_keys = sorted(groups_by_z_key.keys())

    z_levels = []
    cross_section_buckets = []

    index = 0
    while index < len(z_keys):
        z_key = z_keys[index]
        section_labels = groups_by_z_key[z_key]

        z_values = []
        section_index = 0
        while section_index < len(section_labels):
            z_values.append(reference_nodes[section_labels[section_index]][2])
            section_index += 1
        z_level = mean(z_values)

        theta_label_pairs = []
        section_index = 0
        while section_index < len(section_labels):
            node_label = section_labels[section_index]
            x_ref, y_ref, z_ref = reference_nodes[node_label]
            theta = math.atan2(y_ref, x_ref)
            theta_label_pairs.append((theta, node_label))
            section_index += 1

        theta_label_pairs.sort()

        ordered_section_labels = []
        section_index = 0
        while section_index < len(theta_label_pairs):
            ordered_section_labels.append(theta_label_pairs[section_index][1])
            section_index += 1

        z_levels.append(z_level)
        cross_section_buckets.append(ordered_section_labels)
        index += 1

    if len(z_levels) < 2:
        raise RuntimeError('Need at least two z cross-sections in '
                           'INNER_SURFACE')

    return z_levels, cross_section_buckets


# ----------------------------------------------------------------------
# ODB node set lookup
# ----------------------------------------------------------------------

def find_nodeset_in_odb(odb, node_set_name):
    node_set_name_upper = node_set_name.upper()

    if node_set_name_upper in odb.rootAssembly.nodeSets:
        return odb.rootAssembly.nodeSets[node_set_name_upper]

    for instance_name in odb.rootAssembly.instances.keys():
        instance = odb.rootAssembly.instances[instance_name]
        if node_set_name_upper in instance.nodeSets:
            return instance.nodeSets[node_set_name_upper]

    for key in odb.rootAssembly.nodeSets.keys():
        if key.upper() == node_set_name_upper:
            return odb.rootAssembly.nodeSets[key]

    for instance_name in odb.rootAssembly.instances.keys():
        instance = odb.rootAssembly.instances[instance_name]
        for key in instance.nodeSets.keys():
            if key.upper() == node_set_name_upper:
                return instance.nodeSets[key]

    raise RuntimeError('Could not find node set %s in ODB' % node_set_name)


# ----------------------------------------------------------------------
# Read frame coordinates from COORD field
# ----------------------------------------------------------------------

def build_label_to_local_index(inner_surface_labels):
    label_to_local_index = {}
    index = 0
    while index < len(inner_surface_labels):
        label_to_local_index[inner_surface_labels[index]] = index
        index += 1
    return label_to_local_index



# COORD is read at each output frame so the volume is computed from
# the deformed lumen geometry, not the reference mesh.
def frame_coords_from_field(frame, region_nodeset, label_to_local_index,
                            number_of_nodes):
    if 'COORD' not in frame.fieldOutputs:
        raise RuntimeError('COORD field output not found in frame')

    coord_subset = frame.fieldOutputs['COORD'].getSubset(
        region=region_nodeset
    )
    coord_values = coord_subset.values

    current_points = [None] * number_of_nodes

    index = 0
    while index < len(coord_values):
        value = coord_values[index]
        node_label = value.nodeLabel
        if node_label in label_to_local_index:
            local_index = label_to_local_index[node_label]
            current_points[local_index] = (
                value.data[0], value.data[1], value.data[2]
            )
        index += 1

    index = 0
    while index < number_of_nodes:
        if current_points[index] is None:
            raise RuntimeError('Missing COORD for one or more '
                               'INNER_SURFACE nodes in frame')
        index += 1

    return current_points


# ----------------------------------------------------------------------
# Area from one cross-section
# ----------------------------------------------------------------------


# Each cross-section is treated as a polygon in the xy plane and its
# area is computed with the shoelace formula.
def area_from_ordered_section(section_labels, label_to_local_index,
                              current_points):
    # Shoelace polygon area in the xy plane.
    number_of_section_nodes = len(section_labels)
    if number_of_section_nodes < 2:
        return 0.0

    twice_section_area = 0.0
    index = 0
    while index < number_of_section_nodes:
        node_label_1 = section_labels[index]
        node_label_2 = section_labels[(index + 1) % number_of_section_nodes]

        point_1 = current_points[label_to_local_index[node_label_1]]
        point_2 = current_points[label_to_local_index[node_label_2]]

        x_1 = point_1[0]
        y_1 = point_1[1]
        x_2 = point_2[0]
        y_2 = point_2[1]

        twice_section_area += x_1 * y_2 - y_1 * x_2
        index += 1

    section_area = 0.5 * abs(twice_section_area)
    return section_area


# ----------------------------------------------------------------------
# Volume from section areas
# ----------------------------------------------------------------------


# The final lumen volume is the trapezoidal integral of area along z.
def volume_from_sections(z_levels, cross_section_buckets,
                         label_to_local_index, current_points):
    section_areas = []

    index = 0
    while index < len(cross_section_buckets):
        section_area = area_from_ordered_section(
            cross_section_buckets[index],
            label_to_local_index,
            current_points
        )
        section_areas.append(section_area)
        index += 1

    volume = 0.0
    index = 0
    while index < len(z_levels) - 1:
        dz = z_levels[index + 1] - z_levels[index]
        volume += 0.5 * (section_areas[index]
                         + section_areas[index + 1]) * dz
        index += 1

    return volume, section_areas


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

print('Parsing mesh.inp ...')
reference_nodes, inner_surface_labels = parse_inp_nodes_and_nset(
    INP_PATH, NODESET_NAME
)
print('INNER_SURFACE nodes = %d' % len(inner_surface_labels))

print('Building cross-section buckets from reference z ...')
z_levels, cross_section_buckets = build_z_buckets(reference_nodes,
                                                  inner_surface_labels)
print('Cross-sections = %d' % len(z_levels))
print('Selected z-range from INNER_SURFACE:')
print('  z_start = %.6g' % z_levels[0])
print('  z_end   = %.6g' % z_levels[-1])
print('  length  = %.6g' % (z_levels[-1] - z_levels[0]))

index = 0
while index < len(z_levels):
    print('  section %d   z = %.6g   nodes = %d'
          % (index, z_levels[index], len(cross_section_buckets[index])))
    index += 1

label_to_local_index = build_label_to_local_index(inner_surface_labels)

print('Opening ODB ...')
odb = openOdb(ODB_PATH, readOnly=True)
step = odb.steps[STEP_NAME]
nodeset = find_nodeset_in_odb(odb, NODESET_NAME)

number_of_frames = len(step.frames)
print('Frames = %d' % number_of_frames)

time_values = []
volume_values = []
report_stride = max(1, number_of_frames // 10)

print('Computing section areas and volumes ...')
frame_index = 0
while frame_index < number_of_frames:
    frame = step.frames[frame_index]
    current_points = frame_coords_from_field(
        frame,
        nodeset,
        label_to_local_index,
        len(inner_surface_labels)
    )
    volume, section_areas = volume_from_sections(
        z_levels,
        cross_section_buckets,
        label_to_local_index,
        current_points
    )

    time_values.append(frame.frameValue)
    volume_values.append(volume)

    if frame_index % report_stride == 0 or frame_index == number_of_frames - 1:
        print('  frame %4d / %d   t = %.6g   V = %.6g'
              % (frame_index + 1, number_of_frames,
                 frame.frameValue, volume))

    frame_index += 1

odb.close()

V_0 = volume_values[0]
relative_volume_values = []

index = 0
while index < len(volume_values):
    if V_0 != 0.0:
        relative_volume_values.append((volume_values[index] - V_0) / V_0)
    else:
        relative_volume_values.append(0.0)
    index += 1

with open(CSV_OUT, 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['time', 'volume', 'v_rel'])
    index = 0
    while index < len(time_values):
        writer.writerow([time_values[index],
                         volume_values[index],
                         relative_volume_values[index]])
        index += 1

print('--------------------------------------')
print('V0    = %.6g' % V_0)
print('Vmax  = %.6g' % max(volume_values))
print('Vmin  = %.6g' % min(volume_values))
print('Wrote %s' % CSV_OUT)
print('--------------------------------------')
