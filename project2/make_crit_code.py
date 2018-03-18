def modify_template(filename, size_x, size_y, outname):
    with open(filename, 'r') as in_file:
        lines = in_file.readlines()
        for i, line in enumerate(lines):
            if "<<size_x>>" in line:
                lines[i] = line.replace("<<size_x>>", size_x)
            if "<<size_y>>" in line:
                lines[i] = line.replace("<<size_y>>", size_y)

        with open(outname, 'w') as out_file:
            for line in lines:
                out_file.write(line)

def main():
    for L in [30, 60, 120, 240, 360]:
        modify_template('template.cpp', str(L), str(L), 'crit_L_{0}.cpp'.format(L))

if __name__ == '__main__':
    main()
