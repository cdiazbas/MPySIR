import subprocess
import os
import shutil

def main():
    """
    Main function to download, compile, and install the SIR executable.
    Steps:
    1. Remove the existing SIR repository if it exists.
    2. Clone the SIR repository from GitHub.
    3. Change the working directory to the cloned repository.
    4. Compile the SIR code using the provided Python script.
    5. Move the compiled executable to the target directory.
    6. Remove the cloned SIR repository folder.
    Exceptions:
    - Prints an error message if any subprocess command fails or if there are issues with file operations.
    """
    url = 'https://github.com/cdiazbas/SIRcode'
    repo_name = 'SIRcode'
    target_executable_path = '../invDefault/sir.x'

    try:
        # Remove the existing repository if it exists
        if os.path.exists(repo_name):
            shutil.rmtree(repo_name)
            print(f"Removed existing repository folder {repo_name}")

        # Clone the repository
        subprocess.run(['git', 'clone', url], check=True)
        print(f"Successfully cloned repository from {url}")

        # Change to the repository directory
        os.chdir(repo_name)
        print(f"Changed directory to {repo_name}")

        # Compile the code
        with open('compilation_warnings.log', 'w') as log_file:
            subprocess.run(['python', 'createExe.py'], stderr=log_file, check=True)
        print("Successfully compiled the SIR code. Check compilation_warnings.log for details.")

        # Move the executable to the target directory
        os.rename('sir.x', target_executable_path)
        print(f"Moved the executable to {target_executable_path}")

        # Change back to the parent directory and remove the repository folder
        os.chdir('..')
        shutil.rmtree(repo_name)
        print(f"Removed the repository folder {repo_name}")

    except subprocess.CalledProcessError as e:
        print(f"Subprocess error: {e}")
    except FileNotFoundError as e:
        print(f"File operation error: {e}")

if __name__ == "__main__":
    main()
