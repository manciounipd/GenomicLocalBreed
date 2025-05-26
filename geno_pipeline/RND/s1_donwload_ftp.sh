#!/bin/bash


cd Scarico*

# Function to display files in remote directory
display_remote_files() {
    local remote_dir="$1"
    
    # Print current directory
    echo "Files in directory: $remote_dir"
    
    # Get list of files and directories in the remote directory
    files="$(curl --user "$ftp_username:$ftp_password" "ftp://${ftp_server}:${ftp_port}${remote_dir}" -l)"
    
    # Iterate over files and directories
    while IFS= read -r item; do
        # Skip parent directory (.) and current directory (..)
        if [ "$item" != "." ] && [ "$item" != ".." ]; then
            echo "$item"
        fi
    done <<< "$files"
}

# Function to download files from remote directory
download_files() {
    local remote_dir="$1"
    local local_dir="$2"
    
    # Get list of files and directories in the remote directory
    files="$(curl --user "$ftp_username:$ftp_password" "ftp://${ftp_server}:${ftp_port}${remote_dir}" -l)"
    
    # Iterate over files
    while IFS= read -r item; do
        # Skip parent directory (.) and current directory (..)
        if [ "$item" != "." ] && [ "$item" != ".." ]; then
            # Download file
            echo "Downloading file: $item"
            curl --user "$ftp_username:$ftp_password" "ftp://${ftp_server}:${ftp_port}${remote_dir}/${item}" -O
        fi
    done <<< "$files"
}

# FTP server details
ftp_server="fs01.anafi.it"
ftp_port="21"                     
ftp_username="anare"
ftp_password="115oc600"
local_directory="."
remote_directory="/GENOMICA/"

echo "mostra le directory per vefere se funziona..."
display_remote_files "$remote_directory"

echo "controlla quanti file mancano"

ls -d * > old.txt
display_remote_files "/GENOMICA/" | grep "CA"  | sed -e '1d' > new.txt #nb

comm -23  <(sort -k 1 new.txt)  <(sort -k 1 old.txt) > mancano.txt

cat mancano.txt

echo "cominica a scaricare.."

while IFS= read -r i; do
    remote_directory="/GENOMICA/${i}/"
    echo "Creating directory: $i"
    mkdir "$i"
    cd "$i" || exit
    download_files "$remote_directory" .
    cd ..
done < mancano.txt

rm mancano.txt new.txt old.txt

cd ..

echo "renome the folder.."

data=$(date +'%d-%m-%Y')
mv Scarico* "Scarico_$data"