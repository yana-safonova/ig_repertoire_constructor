#!/usr/bin/env python2

import Tkinter, Tkconstants, tkFileDialog
import tkMessageBox


class FileOpenDialog:
    def __init__(self, root, label="LABEL"):
        self.label = Tkinter.Label(root, text=label)
        self.var = Tkinter.StringVar(root)
        self.textbox = Tkinter.Entry(root, textvariable=self.var)

        def browse():
            filename = self.browse()
            if filename:
                self.var.set(filename)

        self.button = Tkinter.Button(root, text="Browse...", command=browse)

        self.label.pack()
        self.textbox.pack()
        self.button.pack()

    def browse(self):
        return tkFileDialog.askopenfilename(initialdir=".",
                                            initialfile=self.var.get(),
                                            title="Select read file",
                                            filetypes=(("FASTQ files", "*.fq"),
                                                       ("FASTA files", "*.fa"),
                                                       ("All files", "*.*")))
    def get(self):
        return self.var.get()


class DirSaveDialog(FileOpenDialog):
    def browse(self):
        return tkFileDialog.askdirectory(initialdir=".",
                                         title="Select output directory",
                                         mustexist=False)


class DropdownMenu:
    def __init__(self, root, options, default_option=0):
        self.var = Tkinter.StringVar(root)
        self.var.set(options[default_option])
        self.option = Tkinter.OptionMenu(root, self.var, *options)

    def get(self):
        return self.var.get()

    def pack(self, *args, **kwargs):
        self.option.pack(*args, **kwargs)


class ThreadMenu:
    def __init__(self, root, default=4):
        self.label = Tkinter.Label(root, text="Threads:")
        self.var = Tkinter.StringVar(root, str(default))
        self.textbox = Tkinter.Entry(root, textvariable=self.var)

    def get(self):
        return int(self.var.get())

    def pack(self):
        self.label.pack()
        self.textbox.pack()


if __name__ == "__main__":
    root = Tkinter.Tk()

    fod = FileOpenDialog(root, label="Merged reads")
    dsd = DirSaveDialog(root, label="Output directory")

    organism = DropdownMenu(root, ["human", "mouse", "pig", "rabbit", "rat", "rhesus_monkey"])
    organism.pack()

    loci = DropdownMenu(root, ["IGH", "IGK", "IGL", "IG", "TRA", "TRB", "TRG", "TRD", "TR", "all"])
    loci.pack()

    threads = ThreadMenu(root, 4)
    threads.pack()

    def run_igrec():
        import os
        home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + "/"
        call = home_directory + "/igrec.py -s %s -o %s -l %s --organism %s -t %d" % (fod.get(), dsd.get(), loci.get(), organism.get(), threads.get())
        return os.system(call)

    run_button = Tkinter.Button(root, text="RUN", command=run_igrec)
    run_button.pack()

    root.mainloop()
