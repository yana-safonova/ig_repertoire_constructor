import glob

def main():
    for filename in glob.glob('../**/*.hpp'):
        with open(filename, "r") as openedFile:
            print filename
            content = openedFile.read()
            changed = False
            if not content.__contains__("#pragma once"):
                content = "#pragma once\n\n" + content
                changed = True
                print "changed"
        if changed:
            with open(filename, "w") as openedFile:
                openedFile.write(content)

if __name__ == '__main__':
    main()
