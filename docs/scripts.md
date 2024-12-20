# Miscellaneous Scripts

## Download and extract tarfile from Google Drive (Python)

```
import os
import tarfile
import urllib.request

url = 'insert_url_here'
cwd = os.getcwd()
dest = os.path.join(cwd, "insert_filename_here")

#create download and extract function
def download_extract():
    try:
        #download from google drive link
        urllib.request.urlretrieve(url, dest)
        print(f'File downloaded successfully: {os.path.abspath(dest)}')

        #open tarfile
        with tarfile.open(dest, 'r:gz') as tar_ref:
            #ensure destination directory exists
            os.makedirs(cwd, exist_ok=True)
            #extract to directory
            tar_ref.extractall(path=cwd)
        print(f'File extracted to: {os.path.abspath(cwd)},spm12')

    except Exception as e:
        print(f'Failed to download and extract file: {e}')
    
    return os.path.join(os.path.abspath(cwd),'spm12')

#call function
download_extract()

```

## Takes a cluster.nii.gz file and outputs borders (contours) around each cluster using fslmaths.

```
#!/bin/bash

#This script takes a cluster.nii.gz file and outputs borders (contours) around each cluster using fslmaths.

cwd=$(pwd)
cluster_nifti="$1" #input filename

#upscale resolution using nearest neighbors interpolation (maintains original cluster values).
flirt -in $cluster_nifti -ref $cluster_nifti -out high_res_clusters.nii.gz -applyisoxfm 0.25 -interp nearestneighbour

#create dilated clusters
fslmaths high_res_clusters.nii.gz -dilF dilated_clusters.nii.gz

#subtract undilated clusters from dilated clusters to form a border, and output cluster_contours.nii.gz
fslmaths dilated_clusters.nii.gz -sub high_res_clusters.nii.gz cluster_contours.nii.gz

#if intermediate file high_res_clusters.nii.gz exists...
if [ -f "${cwd}/high_res_clusters.nii.gz" ]; then
    #remove file
    rm -f high_res_clusters.nii.gz
else
    :
fi

#if intermediate file dilated_clusters.nii.gz exists...
if [ -f "${cwd}/dilated_clusters.nii.gz" ]; then
    #remove file
    rm -f dilated_clusters.nii.gz
else
    :
fi

```

## Unzip all .zip files in current directory and any subdirectories (Bash)
```
#!/bin/bash

#takes path as input. uses cwd if not provided
DIRECTORY="${1:-$(pwd)}"

#use find to locate all .zip files and unzip them in place (in the same directory)
find "$DIRECTORY" -type f -iname "*.zip" | while read -r zipfile; do
  unzip -o "$zipfile" -d "$(dirname "$zipfile")"
done

```

## Extract all .tar files current directory and any subdirectories (Bash)
```
#!/bin/bash

#take path as input or uses cwd if not provided
DIRECTORY="${1:-.}"

#find all .tar files and extract them in place (in same directory)
find "$DIRECTORY" -type f -iname "*.tar*" -print0 | while IFS= read -r -d '' tarfile; do
  echo "Extracting: $tarfile"
  tar -xf "$tarfile" -C "$(dirname "$tarfile")"
  if [ $? -ne 0 ]; then
    echo "Error extracting $tarfile"
  fi
done

```

## Script for converting to BIDS format (Bash)
```
#!/usr/bin/bash

#base directory path var
BASE_DIR="path_to_base_directory"
DEST_DIR="path_to_destination_directory"

#for loop looping through each subdirectory in base
for SUB_DIR in "$BASE_DIR"/*/; do

        #if item is a subdirectory...
        if [ -d "$SUB_DIR" ]; then
                SUBJ_NUM=${SUB_DIR#${BASE_DIR}/LAB_S}
                SUBJ_NUM=${SUBJ_NUM%/}
                echo $SUB_DIR
                echo $SUBJ_NUM
                BIDS_SUB_DIR="${DEST_DIR}/sub-$SUBJ_NUM"
                mkdir -p $BIDS_SUB_DIR

                #define anat path
                ANAT_DIR="${BIDS_SUB_DIR}/anat"

                mkdir -p $ANAT_DIR

                #if the t1 file exists, then move it to anat
                if [ -f "$SUB_DIR/t1.nii.gz" ]; then
                        #move t1 to anat
                        cp $SUB_DIR/t1.nii.gz $ANAT_DIR
                        #rename t1 to proper format
                        echo $ANAT_DIR

                        mv $ANAT_DIR/t1.nii.gz $ANAT_DIR/sub-${SUBJ_NUM}_T1w.nii.gz
                fi

                if [ -f "$SUB_DIR/t2.nii.gz" ]; then
                        #move t2 to anat
                        cp "$SUB_DIR/t2.nii.gz" "$ANAT_DIR"
                        #rename t2
                        mv "$ANAT_DIR/t2.nii.gz" $ANAT_DIR/sub-${SUBJ_NUM}_T2w.nii.gz
                fi
        fi
done


#loops through only directories, -d “$f” checks whether $f is a directory, if true; then action 
for f in *; do
    if [ -d "$f" ]; then
        # $f is a directory
    fi
done
```

## Script for converting to BIDS format (Python)

```
import os
from pathlib import Path

#define base and destination paths 
base_dir_loc = Path("path_to_base_directory")
dest_dir_loc = Path("path_to_destination_directory")

for sub_dir in base_dir_loc.glob('*'):
    if sub_dir.is_dir():
        #convert subdirectory to string and find the subject number as last 5 indexes of the path string
        subj_num = sub_dir.name[-5:]
        
        #create destination subdirectory path with bids format naming
        dest_sub_dir = dest_dir_loc / f"sub-{subj_num}"
        #create destination sub directory
        os.makedirs(dest_sub_dir, exist_ok=True)
        
        #make anat directory
        anat_loc = dest_sub_dir / "anat"
        os.makedirs(anat_loc, exist_ok=True)
        
        #move t1 file to destination anat directory and change name to BIDS format
        t1_files = list(sub_dir.glob("t1.nii.gz")) #find t1 files
        for t1_file in t1_files:
            new_t1_loc = anat_loc/f"sub-{subj_num}_T1W.nii.gz"
            if t1_file.is_file():
                #move and rename file
                t1_file.rename(new_t1_loc)
        t2_files = list(sub_dir.glob("t2.nii.gz")) #find t2 files
        for t2_file in t2_files:
            new_t2_loc = anat_loc/f"sub-{subj_num}_T2W.nii.gz"
            if t2_file.is_file():
                #move and rename file
                t2_file.rename(new_t2_loc) 


```

