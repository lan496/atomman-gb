# atomman-gb

Grain-boundary generator with [Atomman](https://github.com/usnistgov/atomman)

## Gallery

![](examples/gb_13_67.png)

## Usage

```python
from pathlib import Path

import numpy as np
from atomman import Atoms, Box, System

from atomman_gb import CubicGBGenerator, CubicGBInfo


def get_system_fcc() -> System:
    vects = 4.05 * np.eye(3)
    box = Box(vects=vects)
    atype = 1
    pos = [
        [0, 0, 0],
        [0, 0.5, 0.5],
        [0.5, 0, 0.5],
        [0.5, 0.5, 0],
    ]
    atoms = Atoms(atype=atype, pos=pos)

    symbols = ("Al",)
    system = System(atoms=atoms, box=box, symbols=symbols, scale=True)  # set "scale=True" for fractional coordinates!

    return system


if __name__ == "__main__":
    system = get_system_fcc()

    # STGB <001>/(100)
    uvw = (0, 0, 1)
    hkl = (1, 0, 0)

    uvw_str = "".join(map(str, uvw))
    hkl_str = "".join(map(str, hkl))
    root_dir = Path(__file__).resolve().parent / f"Al_stgb_{uvw_str}_{hkl_str}"
    root_dir.mkdir(exist_ok=True)

    gbinfo = CubicGBInfo(uvw=uvw, max_sigma=32)
    for data in gbinfo.datum:
        for planes in data["symmetric_tilt"]:
            plane_mid = planes["plane_mid"]
            if not np.allclose(plane_mid, hkl):
                continue

            gbgen = CubicGBGenerator.make_tilt(
                system, data["rotation"], uvw, data["csl"], planes["plane_1"]
            )
            gb = gbgen.generate(dmin=0.0)[0]  # Fix microscopic DoF for simplicity

            with open(
                root_dir / f"gb_{data['sigma']}_{int(np.degrees(data['theta']))}.atom_dump", "w"
            ) as f:
                atom_dump = gb.dump("atom_dump")
                f.write(atom_dump)
```

## Installation

```shell
git clone git@github.com:lan496/atomman-gb.git
cd atomman-gb
conda create -n atomman-gb python=3.10 pip
conda activate atomman-gb
pip install .
```

## License

atomman-gb is released under a MIT license.

## Development

Installation
```shell
pip install -e ".[dev,docs]" 
```

Document
```shell
sphinx-autobuild --host=0.0.0.0 --port=8000 docs docs_build -v
```

## References

- <https://github.com/wojdyr/gosam>
- <https://github.com/oekosheri/GB_code>
- <https://github.com/spatala/gbpy>
- <https://github.com/ksyang2013/aimsgb>
