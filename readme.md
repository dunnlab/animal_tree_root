# Animal root

A comparison of phylogenetic studies relevant to placing the root of the 
animal phylogeny.

## git LFS

All the `.phy`, `.nex` and `.zip` files in this repo are tracked with [git large file storage](https://git-lfs.github.com/). You'll need to install it following the instructions on the project's website to work with this repo. 

Git LFS keeps the history of the repo cleaner and more performant, while maintaining the data files in the same repo. After installing git LFS and cloning this repo, you can run `git lfs fetch` to get all the lfs-tracked files if they appear to be small plain-text files. Everything else should be automatic.

## scripts directory

To run the python scripts in this directory, you'll need the `animal_root` conda environment. To recreate this environment, install [mini|ana]conda and run:
```
cd scripts
./create_conda_env.sh
```

