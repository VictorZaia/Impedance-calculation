import os

# Define the root directory of your project
root_dir = os.path.dirname(os.path.abspath(__file__))  
output_file = os.path.join(root_dir, "all_code.txt")


# Open the output file to write all contents
with open(output_file, "w", encoding="utf-8") as out_file:
    for foldername, subfolders, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith(".py"):  # Only include Python files
                file_path = os.path.join(foldername, filename)
                out_file.write(f"\n\n# {'='*20} {filename} {'='*20}\n\n")
                
                with open(file_path, "r", encoding="utf-8") as f:
                    out_file.write(f.read())

print(f"All Python files have been copied to {output_file}")
