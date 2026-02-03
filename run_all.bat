@echo off
REM ============================================================================
REM pLIN — Complete Pipeline (Windows)
REM Runs: pLIN assignment → Integration → Figures
REM NOTE: AMRFinderPlus requires WSL or Linux. Use run_all.sh in WSL.
REM ============================================================================

echo ============================================================
echo   pLIN — Complete Analysis Pipeline (Windows)
echo ============================================================
echo.

REM ── Activate virtual environment ─────────────────────────────────────────
if exist .venv\Scripts\activate.bat (
    call .venv\Scripts\activate.bat
    echo   Virtual environment activated.
) else (
    echo   WARNING: .venv not found. Run setup.bat first.
)

REM ── Step 1: Assign pLIN codes ───────────────────────────────────────────
echo.
echo ============================================================
echo   STEP 1/3: Assigning pLIN codes
echo ============================================================
python assign_pLIN.py
echo.

REM ── Step 2: Integrate pLIN + AMR ────────────────────────────────────────
echo ============================================================
echo   STEP 2/3: Integrating pLIN with AMRFinderPlus
echo ============================================================
if exist output\amrfinder\amrfinder_all_plasmids.tsv (
    python integrate_pLIN_AMR.py
) else (
    echo   Skipping: AMRFinderPlus output not found.
    echo   AMRFinderPlus requires Linux/macOS or WSL on Windows.
    echo   Run in WSL: bash run_amrfinder_all.sh
)
echo.

REM ── Step 3: Generate figures ────────────────────────────────────────────
echo ============================================================
echo   STEP 3/3: Generating publication figures
echo ============================================================
if exist output\integrated\pLIN_AMR_integrated.tsv (
    python generate_figures.py
) else (
    echo   Skipping: integrated data not found.
)

echo.
echo ============================================================
echo   Pipeline Complete!
echo ============================================================
echo.
echo   Output files:
echo     output\pLIN_assignments.tsv           - pLIN codes
echo     output\amrfinder\                     - AMR detections
echo     output\integrated\                    - Combined tables
echo     output\figures\                       - Publication figures
echo.
echo ============================================================
pause
