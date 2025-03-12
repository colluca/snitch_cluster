# ----------------------------------------------------------------------
#
# File: gather_data.sh
#
# Last edited: 11.03.2025
#
# Copyright (C) 2025, ETH Zurich and University of Bologna. All rights reserved.
#
# Author: Lorenzo Leone, ETH Zurich
#
# All right reserved
#
# ----------------------------------------------------------------------

#!/usr/bin/env bash
copy_folder() {
    # Check for exactly 3 arguments
    if [ $# -ne 3 ]; then
        echo "Error: You must provide exactly three arguments."
        echo "Usage: $0 <server_name> <scratch_name> <destination_folder_name>"
        echo "---------------------------------------------"
        return 1
    fi

    # Assign arguments to variables
    server_name="$1"
    scratch_name="$2"
    destination_folder_name="$3"

    echo "Server name: $server_name"
    echo "Scratch name: $scratch_name"
    echo "Destination folder name: $destination_folder_name"

    # Create local destination directory if it doesn't exist
    mkdir -p "${destination_folder_name}"

    # Check if the server is reachable via SSH.
    # Check exist code $0 to see if the server exists or is reachable
    ssh -o ConnectTimeout=10 "${server_name}" "exit" &> /dev/null
    if [ $? -ne 0 ]; then
        echo "Error: Server '${server_name}' is not reachable or doesn't exist."
        echo "---------------------------------------------"
        return 1
    fi

    echo "Copying from: /usr/${scratch_name}/${server_name}/colluca/workspace/ISLPED/snitch_cluster..."
    echo "To local folder: ./${destination_folder_name}"

    # Copy synthesis results under the area directory
    cp -r "/usr/${scratch_name}/${server_name}/colluca/workspace/ISLPED/snitch_cluster/nonfree/gf12/fusion/runs/0/reports/15" "./${destination_folder_name}/fusion/"
    # Copy power results under the power directory
    cp -r "/usr/${scratch_name}/${server_name}/colluca/workspace/ISLPED/snitch_cluster/target/snitch_cluster/experiments/frep/power" "./${destination_folder_name}/"

    if [ $? -eq 0 ]; then
      echo "Copy operation completed successfully."
    else
        echo "Error during scp copy operation."
        echo "---------------------------------------------"
      return 1
    fi

    echo "---------------------------------------------"
}



# --------------------------------------------------------
# Call the function for all the experiments to collect
# --------------------------------------------------------

copy_folder larain7  scratch2 zon148dobu
copy_folder larain5  scratch  zonl64dobu
copy_folder larain4  scratch  zonl64fc
copy_folder larain3  scratch2 zonl32fc
copy_folder larain10 scratch  base32fc
