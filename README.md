🔗 Needleman-Wunsch Global Sequence Aligner

A Streamlit web app that performs global sequence alignment using the Needleman-Wunsch dynamic programming algorithm.

What it does

Takes two nucleotide or protein sequences as input
Lets you customise the match reward, mismatch penalty, and gap penalty
Outputs the optimal global alignment with colour-coded matches, mismatches, and gaps
Displays alignment score, match count, mismatch count, gap count, and identity %
Example inputs

Sequence 1	Sequence 2
GATTACA	GCATGCU
ATCGTA	ATCTA
MVLSPADKTNVK	MVHLTPEEKSAVT
Scoring defaults

Parameter	Default
Match reward	+1
Mismatch penalty	-1
Gap penalty	-2
These can be changed in the app.

Run locally

pip install streamlit
streamlit run app.py
Then open http://localhost:8501 in your browser.

Deploy on Streamlit Community Cloud

Push app.py and requirements.txt to a GitHub repo
Go to share.streamlit.io and sign in with GitHub
Click New app → select your repo → set main file to app.py → Deploy
Files

File	Description
app.py	Main Streamlit application
requirements.txt	Python dependencies
README.md	This file
Requirements

Python 3.8+
streamlit >= 1.32.0
No other external libraries needed
