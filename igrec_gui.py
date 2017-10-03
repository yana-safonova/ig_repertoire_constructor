#!/usr/bin/env python2

import Tkinter
import tkFileDialog
from ttkthemes import themed_tk as tk   # Also imports the normal tk definitions, such as Button, Label, etc.
import ttk
import os


home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + "/"


class CompoundControl(object):
    def __init__(self, root):
        self.frame = ttk.Frame(root, height=2)

    def pack(self, *args, **kwargs):
        self.frame.pack(*args, **kwargs)

    def grid(self, *args, **kwargs):
        self.frame.grid(*args, **kwargs)


class FileOpenDialog(CompoundControl):
    def __init__(self, root, label="LABEL"):
        from Tkinter import X, SUNKEN, LEFT, RIGHT, E, W
        CompoundControl.__init__(self, root)

        self.label = ttk.Label(self.frame, text=label)
        self.var = Tkinter.StringVar(self.frame)
        self.textbox = ttk.Entry(self.frame, textvariable=self.var, width=50)

        def browse():
            filename = self.browse()
            if filename:
                self.var.set(filename)

        self.button = ttk.Button(self.frame, text="Browse...", command=browse)

        self.label.pack(side=LEFT)
        self.button.pack(side=RIGHT)
        ttk.Label(self.frame, text=" ").pack(side=RIGHT)
        self.textbox.pack(side=RIGHT)

    def browse(self):
        filetypes = []
        for suffix in ["q", "a"]:
            for compression in ["", ".gz", "bz2"]:
                for prefix in ["f", "fast"]:
                    filetypes.append(("FAST%s files" % suffix.upper(), "*.%s%s%s" % (prefix, suffix, compression)))


        filetypes.append(("All files", "*.*"))
        return tkFileDialog.askopenfilename(initialdir=".",
                                            initialfile=self.var.get(),
                                            title="Select read file",
                                            filetypes=filetypes)

    def get(self):
        return self.var.get()

    def disable(self):
        from Tkinter import DISABLED, NORMAL
        self.textbox.config(state=DISABLED)
        self.button.config(state=DISABLED)

    def enable(self):
        from Tkinter import DISABLED, NORMAL
        self.textbox.config(state=NORMAL)
        self.button.config(state=NORMAL)


class DirSaveDialog(FileOpenDialog):
    def browse(self):
        return tkFileDialog.askdirectory(initialdir=".",
                                         title="Select output directory",
                                         mustexist=False)


class DropdownMenu:
    def __init__(self, root, label, options, default=0):
        from Tkinter import X, SUNKEN, LEFT
        self.frame = ttk.Frame(root, height=2)
        self.label = ttk.Label(self.frame, text=label)
        self.var = Tkinter.StringVar(self.frame)
        if type(default) == int:
            default = options[default]
        self.var.set(default)
        self.option = ttk.OptionMenu(self.frame, self.var, *options)
        self.label.pack(side=LEFT)
        self.option.pack(side=LEFT)

    def get(self):
        return self.var.get()

    def pack(self, *args, **kwargs):
        self.frame.pack(*args, **kwargs)

    def grid(self, *args, **kwargs):
        self.frame.grid(*args, **kwargs)


class ThreadMenu(CompoundControl):
    def __init__(self, root, default=None):
        from Tkinter import X, SUNKEN, LEFT
        CompoundControl.__init__(self, root)

        if default is None:
            import multiprocessing
            cpus = multiprocessing.cpu_count()
            default = cpus - 1 if cpus > 1 else 1

        self.label = ttk.Label(self.frame, text="threads: ")
        self.var = Tkinter.StringVar(self.frame, str(default))
        self.textbox = ttk.Entry(self.frame, textvariable=self.var, width=5)
        self.label.pack(side=LEFT)
        self.textbox.pack(side=LEFT)

        def validate(value):
            import re
            pattern = re.compile("^[0-9]+$")
            return pattern.match(value)

        self.textbox.validatecommand=validate
        self.textbox.validate="all"

    def get(self):
        return int(self.var.get())


class StatusBar(ttk.Frame):
    def __init__(self, master, text=""):
        from Tkinter import SUNKEN, BOTTOM, W, X
        ttk.Frame.__init__(self, master)
        self.label = Tkinter.Label(self, relief=SUNKEN, anchor=W, bg="white")
        self.label.pack(fill=X)
        self.pack(side=BOTTOM, fill=X)
        self.set(text)

    def set(self, text):
        self.label.config(text=text)
        self.label.update_idletasks()

    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()

    def bg(self, color):
        self.label.config(bg=color)


def seticon(root, path):
    import Tkinter
    img = Tkinter.Image("photo", file=path)
    root.tk.call("wm", "iconphoto", root._w, img)


if __name__ == "__main__":
    root = tk.ThemedTk()              # Creates an object for the ThemedTk wrapper for the normal Tk class
    root.set_theme("radiance")
    root.title("Ig Repertoire Constructor")
    root.resizable(False, False)
    root.configure(background="#F6F4F2")

    iconpath = home_directory + "/src/extra/BaseSpace/logos/IgReC_100.png"
    seticon(root, iconpath)

    statusbar = StatusBar(root, "Setup parameters and press RUN")

    inputframe = ttk.LabelFrame(root, text="Input", padding=3)
    inputframe.pack(fill="both", expand="yes", padx=10, pady=5)

    fod = FileOpenDialog(inputframe, label=" ")
    fod_s1 = FileOpenDialog(inputframe, label="Left: ")
    fod_s2 = FileOpenDialog(inputframe, label="Right: ")

    read_type = Tkinter.StringVar()

    from Tkinter import W, E, LEFT, RIGHT, X, TOP
    merged_radio = ttk.Radiobutton(inputframe, text="Merged reads: ",
                                   variable=read_type, value="merged")
    paired_radio = ttk.Radiobutton(inputframe, text="Paired-end reads",
                                   variable=read_type, value="paired-end")

    def reads_select(somevar1, somevar2, mode):
        val = read_type.get()
        if val == "merged":
            fod.enable()
            fod_s1.disable()
            fod_s2.disable()
        elif val == "paired-end":
            fod.disable()
            fod_s1.enable()
            fod_s2.enable()

    read_type.trace("w", reads_select)
    read_type.set("merged")

    barcoded = Tkinter.IntVar(0)
    barcoded_checkbox = ttk.Checkbutton(inputframe, text="barcoded data", variable=barcoded)

    merged_radio.grid(column=0, row=0, stick=W)
    fod.grid(column=1, row=0, stick=E)
    paired_radio.grid(column=0, row=2, stick=W)
    fod_s1.grid(column=0, row=3, stick=E, columnspan=2)
    fod_s2.grid(column=0, row=4, stick=E, columnspan=2)
    barcoded_checkbox.grid(column=0, row=5, stick=W)



    outputframe = ttk.LabelFrame(root, text="Output", padding=3)
    dsd = DirSaveDialog(outputframe, label="Output directory: ")
    dsd.pack(fill=X)
    outputframe.pack(fill="both", expand="yes", padx=10, pady=5)


    paramframe = ttk.LabelFrame(root, text="Parameters", padding=3)
    organism = DropdownMenu(paramframe, label="organism:", options=["human", "mouse", "pig", "rabbit", "rat", "rhesus_monkey"])
    organism.grid(row=0, column=0, padx=0)

    loci = DropdownMenu(paramframe, label="loci:", options=["IGH", "IGK", "IGL", "IG", "TRA", "TRB", "TRG", "TRD", "TR", "all"], default=-1)
    loci.grid(row=0, column=1, padx=20)

    threads = ThreadMenu(paramframe)
    threads.grid(row=0, column=2, padx=0)

    paramframe.pack(fill="both", expand="yes", padx=10, pady=5)


    def input_check():
        if read_type.get() == "merged":
            return fod.get().strip() != ""
        elif read_type.get() == "paired-end":
            return fod_s1.get().strip() != "" and fod_s2.get().strip() != ""

    def output_check():
        return dsd.get().strip() != ""


    def run_igrec():
        if not input_check():
            statusbar.set("Specify input file(s)!")
            statusbar.bg("red")
            return
        if not output_check():
            statusbar.set("Specify output dir!")
            statusbar.bg("red")
            return

        tool = "barcoded_igrec.py" if barcoded.get() else "igrec.py"
        if read_type.get() == "merged":
            data = "-s %s" % fod.get()
        elif read_type.get() == "paired-end":
            data = "-1 %s -2 %s" % (fod_s1.get(), fod_s2.get())

        call = home_directory + "/%s %s -o %s -l %s --organism %s -t %d" % (tool, data, dsd.get(), loci.get(), organism.get(), threads.get())

        statusbar.set("IgReC running...")

        from subprocess import Popen, PIPE, STDOUT
        import shlex
        p = Popen(shlex.split(call), stdin=PIPE, stdout=PIPE, stderr=STDOUT)
        stdoutdata, _ = p.communicate()
        ret_code = p.returncode

        if ret_code == 0:
            statusbar.set("Success")
            statusbar.bg("Green")
        else:
            error = ""
            for line in reversed(stdoutdata.split("\n")):
                if line.strip():
                    error = line.strip()
                    print line
                    break
            statusbar.set("Run failed with code %d: Message: %s" % (ret_code, error))
            statusbar.bg("red")

        print stdoutdata

        return ret_code

    run_button = ttk.Button(root, text="RUN", command=run_igrec)
    run_button.pack(padx=10, pady=5)

    root.mainloop()
