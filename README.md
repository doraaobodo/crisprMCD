# crisprMCD CLI
Command line interface (CLI) to run robust minimum covariant determinant (MCD) analyses on large-scale CRISPR studies involving multiple CRISPR knockout screens using R. This project provides a simple cli (for any OS) for an interactive R workflow for MCD analysis.
After cloning or downloading the project folder, run the launcher for your operating system (run.bat for Windows OS Launcher and .sh for MAC/Linus Launcher) to start.

## Requirements

- R must be installed on your machine
- No R IDE is required
- The launcher will try to find `Rscript` automatically

If R is not installed, please install it first:

- CRAN: https://cran.r-project.org/

---

## Important Project Files

```text
main.R        # main interactive R workflow
run.bat       # Windows launcher
run.sh        # macOS/Linux launcher
data/         # example input data
```

---

## Download the Project

### Option 1: Download from GitHub as a ZIP
1. Open the repository page on GitHub
2. Click **Code**
3. Click **Download ZIP**
4. Extract the ZIP file to a location on your computer

### Option 2: Clone with git

```bash
git clone https://github.com/doraaobodo/crisprMCD.git
cd crisprMCD
```

---

## How to Run

### Windows

1. Open the project folder

```bash
cd crisprMCD
```
2. Double-click `run.bat`

The launcher will:
- look for the newest available R installation
- run `main.R`
- report any errors clearly

If double-clicking does not work, open Command Prompt in the project folder and run:

```bat
run.bat
```

---

### macOS / Linux

1. Open a terminal
2. Navigate to the project folder

```bash
cd crisprMCD/
```
3. Make the launcher executable:

```bash
chmod +x run.sh
```
4. Make sure we can access osascript

Go to:

System Settings -> Privacy & Security -> Automation

and ensure:

* Terminal
* iTerm
* R
* RStudio

have permission to control Finder/System Events.

5. Run:

```bash
./run.sh
```

The launcher will:
- look for a usable R installation
- run `main.R`
- report any errors clearly

** Be Aware **

We use the XQuartz (X11) editor to edit columns of input files. 
XQuartz (X11) terminal on macOS often fails to render or accept special characters because it defaults to a limited character set.
To change column header names, type in a cell, copy and paste edit into column header cell. 
Clip below demonstrates:
---

### Alternative

If the launcher scripts do not work, you can run the project directly with:

```bash
Rscript main.R
```

---

---

## Example Data

Example input files are included in the `data/` folder so users can test the workflow.
This is a MAGECK-MLE results file.
Example:

```text
data/xiao_2018_mle.csv
```

Use this file and specify MLE when prompted by the application.

---

## Notes

- The launcher only starts the application
- All interactive behavior is handled in `main.R`
- If required R packages are missing, R may prompt you or stop with an error.

---

## Troubleshooting

### `Rscript` was not found
R is either not installed or not accessible on your machine.

Install R from:
https://cran.r-project.org/

Then try again.

### `main.R` was not found
Make sure you are running the launcher from inside the full downloaded project folder and that `main.R` is present.

### Permission denied on macOS/Linux
Run:

```bash
chmod +x run.sh
```

and then try again.

---

## Contact / Support

If you encounter issues, please open an issue on this GitHub repository and include:
- your operating system
- your R version
- the error message shown in the launcher