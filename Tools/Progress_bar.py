import sys

def progress_bar(current, total, prefix='Progress', length=40, fill='█', print_end="\r"):

    percent = ("{0:.1f}").format(100 * (current / float(total)))
    
    # Calculate the number of filled positions
    filled_length = int(length * current // total)
    
    # Create the progress bar
    bar = fill * filled_length + '-' * (length - filled_length)
    
    # Print the progress bar, overwriting the last line
    sys.stdout.write(f'\r{prefix} |{bar}| {percent}%')

    sys.stdout.flush()