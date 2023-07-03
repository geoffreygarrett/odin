import os
import sys

root = 'include/odin/'
output_dir = sys.argv[1]  # Get output directory from command line argument

for directory, _, files in os.walk(root):
    headers = [f for f in files if f.endswith('.hpp')]
    if headers:
        # Write to the output directory specified as a command line argument
        with open(os.path.join(output_dir, f'{os.path.basename(directory)}.hpp'), 'w') as f:
            for header in headers:
                f.write(f'#include "{header}"\n')
