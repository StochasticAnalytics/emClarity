#!/bin/bash

# If a binary is accidently committed, the cleanest removal is this.

git filter-branch --index-filter 'git rm --cached --ignore-unmatch bin/emClarity_*' -- --all
rm -Rf .git/refs/original
rm -Rf .git/logs/
git gc --aggressive --prune=now
