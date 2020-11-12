This directory contains some convenient tools and scripts, mostly for developers.

* apply_astyle.py: apply astyl to all C and C++ files in the current directory and its subfolders.

    Usage:
    
    ```
    cd siconos_source_dir
    python scripts/devel/apply_astyle.py
    ```
    
* update_license_header.py: update all Siconos files with the proper cartdrige describing the license.
  
  Usage:
  
  edit the file to set the string that must be replaced in all files and run (from siconos source dir):
      
  ```
  python scripts/devel/update_license_header.py
  ```

