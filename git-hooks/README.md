# GIT Hooks

This directory contains hook scripts that are useful for this project.

**Please copy any hooks into the `.git/hooks` directory**

This cannot be done automatically - GIT does not allow doing this for good reasons, because hooks can contain any command, so installing hooks automatically would be a severe security risk. Please feel free to theck the content of the pre-commit hook to understand what exactly it is doing.

## `pre-commit`

This script increases the version number of the nirvanaparser library
with each commit.
