#!/usr/bin/python3

import os

file_wav = os.getenv('file_wav')
file_wav_prm = os.getenv('file_wav_prm')
file_image = os.getenv('file_image')
nx = int(os.getenv('nx'))
ny = int(os.getenv('ny'))
electron_energy = float(os.getenv('electron_energy'))
Nlp = int(os.getenv('Nlp'))
Nabr = int(os.getenv('Nabr'))

with open(file_wav_prm, "w") as f:
    print(f"'{file_wav}'", file=f)
    print(f"{nx} {ny}", file=f)
    print(f"{0.1} {0.1}", file=f)
    print(f"{electron_energy}", file=f)
    print(f"{0}", file=f)
    print(f"'{file_image}'", file=f)
    print(f"{nx} {ny}", file=f)
    print(f"{0} {1.0} {1.0} {0.0}", file=f)
    print(f"{0}", file=f)
    print(f"{0.1}", file=f)
    print(f"{0.0} {0.0}", file=f)
    print(f"{0.0}", file=f)
    print(f"{1}", file=f)
    print(f"{0} {1.0}", file=f)
    print(f"{0} {1.0}", file=f)
    print(f"{0} {1.0} '{'mtf.dat'}'", file=f)
    print(f"{0} {1.0} {1.0} {0.0}", file=f)
    print(f"{Nabr}", file=f)
    for i in range(Nabr):
        print(f"{0} {0.0} {0.0}", file=f)
    print(f"{250.0} {0.03}", file=f)
    print(f"{0.0} {0.0}", file=f)
    print(f"{Nlp}", file=f)
    for i in range(Nlp):
        print(f"{1}", file=f)
        print(f"{1}", file=f)
        print(f"{1}", file=f)
        print(f"{0.0} {10.0} {21}", file=f)
        print(f"''", file=f)
