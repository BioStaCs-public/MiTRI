# GitHub upload commands

These commands are for the repository owner. Run them in PowerShell from the repository root.

```powershell
cd "<repo-root>"
```

Moves into the local repository folder.

```powershell
git status
```

Checks whether there are uncommitted local changes.

```powershell
git branch -M main
```

Ensures the branch name is `main`.

```powershell
git remote -v
```

Checks the configured GitHub remote.

```powershell
git remote set-url origin https://github.com/BioStaCs-public/MiTRI.git
```

Sets the target GitHub repository.

```powershell
git push -u origin main
```

Uploads the committed code to GitHub and links the local `main` branch to the remote `main` branch.

If GitHub asks for login, use the GitHub account that owns or can write to `BioStaCs-public/MiTRI`. If password login is rejected, use a GitHub personal access token as the password.
