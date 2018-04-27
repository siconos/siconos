"""Python utility to create sphinx rst file from png generated
by doxygen (class diagrams ...)

- scan doxygen ouput (html) path
- create an image entry for each class_inh*.png file in a rst file

"""
import glob
import os

# Scan doxygen output path and create a list with
# files matching requirements
html_doxygen_output_path = '@DOC_ROOT_DIR@/html/doxygen'
outputfile = '@CMAKE_CURRENT_BINARY_DIR@/sphinx/reference/class_diagrams.rst'
class_diagram_match = 'inherit_graph*.png'

header = '.. _api_class_diagrams:\n\n'
header += 'Siconos API - Classes diagrams\n'
header += '==============================\n\n'

files = glob.glob(os.path.join(html_doxygen_output_path, class_diagram_match))
# we need relative paths for sphinx ...
realfiles = [f.split('@DOC_ROOT_DIR@/html')[1] for f in files]
realfiles = ['/..' + f for f in realfiles]

file = open(outputfile, 'w')
file.writelines(header)
#params = [':height: 190 px', ':class: gallery']
params = [':class: gallery']
img_prefix = '.. image:: '

for f in realfiles :
    line = img_prefix + f
    for p in params:
        line += '\n\t' + p
    line += '\n\n'
    file.writelines(line)

file.close()
