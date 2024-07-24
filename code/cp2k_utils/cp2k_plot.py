"""
2024/07/24

Plot cp2k generated output from simple md runs

"""
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

def read_cp2k_file(file_path):
    # Read the file into a pandas DataFrame, skipping the first row which is a comment
    df = pd.read_csv(file_path, delim_whitespace=True, skiprows=1,
                     names=['Step', 'Time_fs', 'Kin_a.u.', 'Temp_K', 'Pot_a.u.', 'ConsQty_a.u.', 'UsedTime_s'])
    return df

def plot_and_save(df, output_folder):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Plot kinetics energy over time
    plt.figure(figsize=(10, 6))
    plt.plot(df['Time_fs'], df['Kin_a.u.'], marker='o', linestyle='-', color='b')
    plt.title('Kinetics Energy over Time')
    plt.xlabel('Time [fs]')
    plt.ylabel('Kinetics Energy [a.u.]')
    plt.grid(True)
    plt.savefig(os.path.join(output_folder, 'kinetics_energy.png'))
    plt.close()

    # Plot temperature over time
    plt.figure(figsize=(10, 6))
    plt.plot(df['Time_fs'], df['Temp_K'], marker='o', linestyle='-', color='r')
    plt.title('Temperature over Time')
    plt.xlabel('Time [fs]')
    plt.ylabel('Temperature [K]')
    plt.grid(True)
    plt.savefig(os.path.join(output_folder, 'temperature.png'))
    plt.close()

    # Plot potential energy over time
    plt.figure(figsize=(10, 6))
    plt.plot(df['Time_fs'], df['Pot_a.u.'], marker='o', linestyle='-', color='g')
    plt.title('Potential Energy over Time')
    plt.xlabel('Time [fs]')
    plt.ylabel('Potential Energy [a.u.]')
    plt.grid(True)
    plt.savefig(os.path.join(output_folder, 'potential_energy.png'))
    plt.close()


FILE_PATH = sys.argv[1]
OUTPUT_FOLDER = sys.argv[2]

if not os.path.isfile(FILE_PATH):
    print(f'File {FILE_PATH} not found.')
    sys.exit(1)

#check arguments
if len(sys.argv) != 3:
    print(f'Usage: {sys.argv[0]} <file_path> <output_folder>')
    sys.exit(1)

def main():
    # Read data from the file
    df = read_cp2k_file(FILE_PATH)
    
    # Plot and save the graphs
    plot_and_save(df, OUTPUT_FOLDER)
    
    print(f'Plots saved in {OUTPUT_FOLDER} folder.')

if __name__ == '__main__':
    main()

