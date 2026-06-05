# Herwig QCD Instantons

## Installing `MamboPhasespace` into a Herwig source tree

Let `HERWIG_SRC` be the top-level Herwig source directory, for example the directory that contains `configure` and `MatrixElement/`.

1. Copy the phase-space source files:

```sh
cp Phasespace/MamboPhasespace.cc "$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/"
cp Phasespace/MamboPhasespace.h  "$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/"
```

2. Edit `$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/Makefile.am`:

Add `MamboPhasespace.h` to `ALL_H_FILES` and `MamboPhasespace.cc` to `ALL_CC_FILES`.

For an already-configured build tree, either rerun the autotools/configure step, or make the same additions in the generated `Makefile.in` and `Makefile`.

3. Rebuild and install Herwig:

```sh
cd "$HERWIG_SRC"
make -j"$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
make install
```

4. Check that Herwig can instantiate the new class:

```sh
cat >/tmp/mambo-phasespace-smoke.in <<'EOF'
cd /Herwig/MatrixElements/Matchbox/Phasespace
create Herwig::MamboPhasespace MamboPS
EOF

Herwig read /tmp/mambo-phasespace-smoke.in
```

The final command should exit without an error.
