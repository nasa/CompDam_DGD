# Advanced debugging
Using an interactive debugger helps to identify issues in the Fortran code. Abaqus knowledge base article QA00000007986 describes the details involved. The following is a quick-start guide for direct application to CompDam.

Several statements for debugging need to be uncommented in the environment file. Follow these steps:
1. Copy your system environment file to your local working directory. For the example below, copy the environment file to the `tests` directory.
2. Edit the local environment file: uncomment lines that end with `# <-- Debugging`, `# <-- Debug symbols`, and `# <-- Optimization Debugging`

Run the job with the `-debug` and `-explicit` arguments. For example:
```
$ abaqus -j test_C3D8R_fiberTension -user ../for/CompDam_DGD.for -double both -debug -explicit
```

This command should open the [Visual Studio debugging software](https://msdn.microsoft.com/en-us/library/sc65sadd.aspx) automatically. Open the source file(s) to debug. At a minimum, open the file with the subroutine entry point `for/CompDam_DGD.for`. Set a break point by clicking in the shaded column on the left edge of the viewer. The break point will halt execution. Press <kbd>F5</kbd> to start the solver. When the break point is reached, a yellow arrow will appear and code execution will pause. Press <kbd>F5</kbd> to continue to the next break point, press <kbd>F11</kbd> to execute the next line of code following execution into function calls (Step Into), or press <kbd>F10</kbd> to execute the next line of code but not follow execution into function calls (Step Over).

To stop execution, close the Visual Studio window. Choose stop debugging and do not save your changes.

[More tips on debugging Fortran programs from Intel](https://software.intel.com/en-us/articles/tips-for-debugging-run-time-failures-in-intel-fortran-applications).

In case you must use remote ssh/scp access to run CompDam, a neat trick is to use the `Commands -> Keep Remote Directory up to Date...` option in [WinScp](https://winscp.net/eng/index.php). This feature can be used during development so that the CompDam files are edited locally and then automatically synced on a remote server for testing. The following mask (`Keep Remote Directory up to Date... -> Transfer Settings... -> File mask:`) can be used for syncing the source code: `*.for; *.py; *.inp; *.props; *.inc; Makefile; kind_map | tests/testOutput/; .git/; .vscode/`.
