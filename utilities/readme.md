# CompDam Utilities

This directory includes utilities that have been developed in parallel with CompDam. The utilities are organized into directories. Installation and usage information are included for each utility.

## Listing of utility tools
1. `ramberg-osgood-fit`: A python script to fit a Ramberg-Osgood curve to stress-strain data.
2. `link-checker`: A standalone python script to check internal and external links in markdown files.
3. `abaqus-python-addpkg`: Removed in CompDam 2.7 since modern versions of Abaqus/CAE include pip. See [How to Install Python Packages in Abaqus](https://tecnodigitalschool.com/how-to-install-python-packages-in-abaqus/).
4. `meshtools`: Removed in CompDam 2.7, may be added back in a future release.

## Link checker usage

Run from the repository root:

```bash
python utilities/link-checker/check_markdown_links.py .
```

The script exits with code `1` when one or more broken/dead links are found, which allows it to be used in CI jobs.
