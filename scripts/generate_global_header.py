import os
import sys


def main():
    output_file = sys.argv[1]
    src_files = sys.argv[2:]

    with open(output_file, 'w') as out_f:
        current_dir = None
        for src_file in src_files:
            directory, filename = os.path.split(src_file)
            # we only consider the directory inside "odin/"
            sub_dir = directory.split("/odin/")[-1]
            if current_dir != sub_dir:
                out_f.write(f"\n// Headers from {sub_dir}\n")
                current_dir = sub_dir
            relative_path = f"odin/{sub_dir}/{filename}"
            out_f.write(f'#include "{relative_path}"\n')


if __name__ == "__main__":
    main()
