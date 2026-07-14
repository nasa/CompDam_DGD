# Checklist for contributions to CompDam

Pull requests (PR) to CompDam are welcome and much appreciated. Here is a checklist to go through before submitting a PR:
1. The new code is formatted consistently with the existing code base
2. One or more tests is included that demonstrate and verify the new functionality
3. The new functionality is documented in the README as appropriate
4. All tests are passing (please attach the test report to the PR)
5. Spell check
6. Run the link-checker in the utilities section: `python utilities/link-checker/check_markdown_links.py .`
7. Generate pes files for examples: `python dcb_runner.py --genPes`
8. Update the section 'Citing CompDam' in the README with the new proposed version number
