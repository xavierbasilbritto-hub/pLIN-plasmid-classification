@echo off
REM ============================================================================
REM pLIN — Complete Pipeline (Windows)
REM Runs: pLIN assignment → Integration
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
echo   STEP 1/2: Assigning pLIN codes
echo ============================================================
python assign_pLIN.py
echo.

REM ── Step 2: Integrate pLIN + AMR ────────────────────────────────────────
echo ============================================================
echo   STEP 2/2: Integrating pLIN with AMRFinderPlus
echo ============================================================
if exist output\amrfinder\amrfinder_all_plasmids.tsv (
    python integrate_pLIN_AMR.py
) else (
    echo   Skipping: AMRFinderPlus output not found.
    echo   AMRFinderPlus requires Linux/macOS or WSL on Windows.
    echo   Run in WSL: bash run_all.sh
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
echo.
echo   Launch the GUI:
echo     streamlit run plin_app.py
echo.
echo ============================================================
pause
