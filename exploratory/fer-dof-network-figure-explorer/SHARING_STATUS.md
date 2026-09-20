# Sharing details

**FER-DOF Network Figure Explorer**  
**Prototype; not yet carefully vetted.** Prepared 20 September 2026.

The package belongs in the existing private [FER paper GitLab repository](https://gitlab.caltech.edu/nolan_lab/feronia-cell-type-specific), on branch `codex/fer-dof-figure-explorer`, under `exploratory/fer-dof-network-figure-explorer`.

- [Start here on GitLab](https://gitlab.caltech.edu/nolan_lab/feronia-cell-type-specific/-/tree/codex/fer-dof-figure-explorer/exploratory/fer-dof-network-figure-explorer)
- [Named PDF with repository source links](FER_DOF_Network_Figure_Explorer.pdf)
- [Source, input and code guide](SOURCES.md)
- [Unresolved validation issues](VALIDATION_STATUS.md)

Chia-Yu Chen (`cchen10`) has existing Owner access inherited from `Nolan_lab`, verified in the repository's Members view on 20 September 2026. Repository visibility was verified as Private. No permissions were changed.

The prototype branch starts from upstream commit `5649a9485de0d50d03c8135c5ff254bdaa8d0a1a`. The original paper analysis and `master` are preserved. Use GitLab's commit view for the exact saved package revision. Preparing this branch does not constitute scientific approval or a change to the manuscript's selected figures.

The share PDF adds a named cover and replaces 32 local-file annotations with links to packaged files on this branch. The original six-page proof remains unchanged. [PDF link checks](provenance/SHARE_PDF_LINK_CHECK.json) record the source and output hashes. Readers need repository access for these links.

To obtain the code and inputs locally:

```bash
git clone --single-branch --branch codex/fer-dof-figure-explorer \
  git@gitlab.caltech.edu:nolan_lab/feronia-cell-type-specific.git
cd feronia-cell-type-specific/exploratory/fer-dof-network-figure-explorer
python3 scripts/verify_snapshot.py
```

The local verification checks saved-file identity and package completeness. It does not validate the GRNs or resolve the source findings. No message has been sent to Chia-Yu on behalf of the user.
