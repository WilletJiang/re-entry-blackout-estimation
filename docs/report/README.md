# Entry Blackout Estimation Report

This directory contains the LaTeX source code for the Entry Blackout Estimation technical report.

## Compiling

To compile the report, use `xelatex` (required for advanced font support):

```bash
cd docs/report
xelatex main.tex
```

## Customization

To customize the font paths or add your own, edit `docs/report/structure.sty`. The font loading logic is contained within `\IfFontExistsTF` blocks to prevent compilation errors on systems without the premium assets.

> **Note on Open Source**: The commercial font files (`.otf`/`.ttf`) are **NOT** included in this repository to comply with licensing agreements. The PDF output may be distributed, but the raw font files must remain on the licensed user's local machine.
