import tkinter as tk
import os
from tkinter import filedialog, messagebox

from main import run_code
#imports run_code function from main.py

def select_file(file_var):
    file_path = filedialog.askopenfilename(parent= root, title="Select a file")
    if file_path:
        file_var.set(file_path) #store file path in file_var (which is parsed into the function)
    else:
        file_var.set("")
#opens file dialog and stores file path in a variable (file_var), else sets to an empty string (i.e. if canceled)

def file_processing():
    mut_path = file_var_gen.get()
    meta_path = file_var_meta.get()
    gtf_path = file_var_ann.get()
    output_path = file_var_out.get()
#gets file paths from GUI and stores in path variables

    files = {"Genetic data": mut_path, "Annotation": gtf_path, "Metadata": meta_path}
    missing_files = [name for name, path in files.items() if not path or not os.path.isfile(path)] #used to check if files exist
    if missing_files:
        messagebox.showwarning("Missing file", "Please select valid files")
        return #stops if files are missing
    try:
        run_code(mut_path = mut_path, gtf_path = gtf_path, meta_path = meta_path, output_path= output_path) #passes file paths to run_code function in main file
        messagebox.showinfo("Process", "Successful")
    except Exception as e:
        messagebox.showerror("Error processing files", str(e)) #shows error message if fails
root = tk.Tk()
root.title("Home")
root.geometry("400x400")
#creates app window, sets gui window size and title

file_var_gen = tk.StringVar()
file_var_meta = tk.StringVar()
file_var_ann = tk.StringVar()
file_var_out = tk.StringVar()
#variables for file paths

#Button genetic data
tk.Button(root, text="Genetic data file path", command = lambda: select_file(file_var_gen)).pack()
tk.Entry(root, textvariable=file_var_gen).pack()
#Button metadata
tk.Button(root, text="Metadata file path", command = lambda: select_file(file_var_meta)).pack()
tk.Entry(root, textvariable=file_var_meta).pack()
#Button annotation
tk.Button(root, text="Annotation file path", command = lambda: select_file(file_var_ann)).pack()
tk.Entry(root, textvariable=file_var_ann).pack()
#Output path
tk.Button(root, text="output file path", command = lambda: select_file(file_var_out)).pack()
tk.Entry(root, textvariable=file_var_out).pack()
#generally, these buttons call select_file with the associated variable (i.e. file_var_gen)
#although, output file path should be updated from a select file to a save file dialog in the next iteration

#Button process files
tk.Button(root, text="Process files", command = file_processing).pack()
#this button runs the file_processing function which is linked to the run_code function in main.py file
root.mainloop()


