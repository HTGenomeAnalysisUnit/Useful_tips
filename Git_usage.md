# GIT useful commands

- [GIT useful commands](#git-useful-commands)
  - [Init a new repo in a folder and connect to gitlab/github](#init-a-new-repo-in-a-folder-and-connect-to-gitlabgithub)
  - [Create a new branch](#create-a-new-branch)
  - [Delete a branch](#delete-a-branch)
    - [locally](#locally)
    - [remotely](#remotely)
  - [Merge a branch into main](#merge-a-branch-into-main)
  - [Change remote URL for an existing repo](#change-remote-url-for-an-existing-repo)
  - [Undo a commit](#undo-a-commit)
    - [Undoing a Specific Commit (That Has Been Pushed)](#undoing-a-specific-commit-that-has-been-pushed)
    - [Undoing Your Last Commit (That Has Not Been Pushed)](#undoing-your-last-commit-that-has-not-been-pushed)
    - [Undoing Local Changes That Have Been Committed (But Not Pushed)](#undoing-local-changes-that-have-been-committed-but-not-pushed)

## Init a new repo in a folder and connect to gitlab/github
```bash
git init
git add --all
git commit -m "first commit"
git branch -M main
git remote add origin $git_repo
git push -u origin main
```

## Create a new branch
```bash
git checkout -b new-branch
git push --set-upstream origin new-branch
```

## Delete a branch
### locally
`git branch -d localBranchName`

### remotely
`git push origin --delete remoteBranchName`

## Merge a branch into main
```bash
git checkout main
git merge new-branch
git push
```

## Change remote URL for an existing repo
`git remote set-url origin https://git-repo/new-repository.git`

## Undo a commit
### Undoing a Specific Commit (That Has Been Pushed)

If you have one specific commit you want to undo, you can revert it as follows:

1. From your git project folder run `git status` and make sure you have a clean working tree.
2. Find the hash for the commit you want to undo (each commit has a unique hash like 2f5451f). Ypu can find commits hashes in 2 ways:
  - In the commit history on the GitHub / GitLab
  - In your terminal using the command `git log --oneline`
3. Once you know the hash for the commit you want to undo, run the following command: 

`git revert [hash-key] --no-edit`
  
  **NOTE:** The --no-edit option prevents git from asking you to enter in a commit message. If you don't add that option, you'll end up in the VIM text editor. To exit VIM, press : to enter command mode, then q for quit, and finally hit Return (Mac) or Enter (Windows).
  This will make a new commit that is the opposite of the existing commit, reverting the file(s) to their previous state as if it was never changed.
4. If working with a remote repo, you can now push those changes: `git push`

### Undoing Your Last Commit (That Has Not Been Pushed)

If you made a mistake on your last commit and have not pushed yet, you can undo it. 
Your changes remain in place, so you can make any additional changes or add any missing files. You can then make a new commit.

From your git project folder run the command `git reset --soft HEAD~`

**TIP:** Add a number to the end to undo multiple commits. For example, to undo the last 2 **unpushed** commits run `git reset --soft HEAD~2`
**NOTE:** `git reset --soft HEAD~` is the same as `git reset --soft HEAD^`


### Undoing Local Changes That Have Been Committed (But Not Pushed)

If you have made local commits that you don't like, and they have not been pushed yet you can reset things back to a previous good commit.
All the files will be reverted to their status before the commit.

1. From your project folder find the hash for the last good commit (the one you want to revert back to) from the GitHub website or using `git log --oneline`
2. Using the hash key you can do
  - `git reset [hash-key]` : commits will be removed, changes will appear as uncommitted but files are left untouched    
  - `git reset --hard [hash-key]`: undo the commits and through away the code, your files are reverted to the status of the commit
