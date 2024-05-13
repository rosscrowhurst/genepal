#!/usr/bin/env bash

perl -p -e 's/^(>\S+).*$/$1/' \
    modules/kherronism/braker3/tests/test.fa
