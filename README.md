# G4ShapePredictor (G4SP)
G4ShapePredictor is an application designed to accurately predict DNA G4 folding topologies in potassium (K+) buffer based on putative quadruplex sequences. For a detailed examination of the source code used in the paper, please see section on Supplementary Code.

Liew, D., Lim, Z. W., & Yong, E. H. (2024). Machine learning-based prediction of DNA G-quadruplex folding topology with G4ShapePredictor. Scientific Reports, 14(1), 24238.

# Installation
## G4ShapePredictor.exe (Windows OS only)
If you're using Windows OS, the simplest way to install G4ShapePredictor is to download the .exe file ([Google Drive](https://drive.google.com/file/d/1Cere2jX_zUgypFrG4FWjMmpy0yBjqASH/view?usp=sharing), [OneDrive](https://entuedu-my.sharepoint.com/:u:/g/personal/liew0207_e_ntu_edu_sg/EWva8coMN1ROvsZfTWumbHABUtiK4jk7lWsIOJM-nsM_KQ?e=k1LcmO), [MediaFire](https://www.mediafire.com/file/if9f4d713b4tiew/G4ShapePredictor.zip/file)). The files are hosted on an external site due to GitHub's upload limit.

1. Download G4ShapePredictor.zip
2. Unzip G4ShapePredictor.zip
3. Right click G4ShapePredictor.exe > Run as Administrator

If there are errors, include the G4ShapePredictor folder as an exception in Windows Defender, firewall, or any antivirus software you might be using.

## Running G4ShapePredictor.py using Python (Windows OS, MacOS, Linux)

### WindowsOS
1. Download and install [Miniconda](https://docs.anaconda.com/miniconda/)
2. Download (or clone) this repository. You may download via: https://github.com/donn-liew/G4ShapePredictor > green "Code" button > download ZIP. Unzip G4ShapePredictor-main.zip.
3. Open Anaconda Prompt from the start menu
4. Run the following command:
```bash
conda env create -f path/to/G4ShapePredictor-main/environment.yml
```
Example: if G4ShapePredictor-main.zip was extracted to my desktop, I would run the following command:
```bash
conda env create -f C:/Desktop/G4ShapePredictor-main/environment.yml
```
5. In Anaconda Prompt, activate the virtual environment by running:
```bash
conda activate github-g4predictor
```
6. Run G4ShapePredictor.py:
```bash
python path/to/G4ShapePredictor-main/g4sp application code/G4ShapePredictor.py
```

### MacOS and Linux
1. Download [Miniconda](https://docs.anaconda.com/miniconda/)
2. Open a terminal window, navigate to the directory containing the downloaded script, and run it with bash Miniconda3-latest-MacOSX-x86_64.sh (adjust the script name as needed for your OS and version). Follow the on-screen prompts.
3. Download (or clone) this repository. You may download via: https://github.com/donn-liew/G4ShapePredictor > green "Code" button > download ZIP.
4. Navigate to the downloaded repository, unzip the file and navigate to the unzipped file:
```bash
cd /path/to/downloaded/repository
unzip G4ShapePredictor-main.zip
cd G4ShapePredictor-main
```
5. Create and activate virtual environment from environment.yml:
```bash
conda env create -f environment.yml
conda activate github-g4predictor
```
6. Run G4ShapePredictor.py:
```bash
python g4sp application code/G4ShapePredictor.py
```
<!-- 
# Supplementary Code
1. Code is written in Python programming language and is found in g4sp supplementary > G4ShapePredictor_Paper_Supplementary.py
2. Packed neatly as a class function G4Data
-->
