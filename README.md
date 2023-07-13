# discrete-log

This repo contains the source code to our submission

> [Dlog is Practically as Hard (or Easy) as DH - Solving Dlogs via DH Oracles on EC Standards](https://eprint.iacr.org/2023/539)

## Files

### `src-codebook`

This folder contains the C code to create the codebook, given a base curve
with appropriate auxiliary curve and generator points.

```
Usage:
  codebook-gen p E_a4 E_a6 Px Py q I_a4 I_a6 Rx Ry order_R factor name
```

### `src-sage`

This folder contains the code to search for auxiliary curves (`curve_finder.sage`)
and to run the implementation with an oracle (`maurer.sage`).

Parameters for base curves and auxiliary curves are given in `curve_list.py`.

Note that the python code has an additional dependency on [`prtpy`](https://pypi.org/project/prtpy/).
Make sure that it is available to sage (e.g. by running `pip install prtpy` or using the [nix flake](https://determinate.systems/posts/nix-run) in this repository).
