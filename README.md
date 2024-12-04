## Systems Architecture Project, FIT ##
> Need scipy <br />
> run this in your cmd: <br />
> python -m pip install scipy <br />


‎

# To start using (How I do it anyways)
- **Download _git bash_**
- **Open _git bash_ and navigate to foulder that you would like this to live**
- **Create a folder in your workspace with the name _'SystemsArchitectureFITProject'_ and go inside that folder**
```
mkdir SystemsArchitectureFITProject && cd SystemsArchitectureFITProject
```
- **Initialize a local git repo**
```
git init -b main
```
- **Grab the _REMOTE-URL_ (found in the 'code' tab above or in address bar) and connect your local repo to this one**
```
git remote add origin REMOTE-URL
```
- **Optional:** You can verifiy that the repo is added correctly by checking the output of this command matches
```
git remote -v
```
- **Lastly, do a _git pull_ to update your folder with the files on the branch**
```
git pull origin main
```
Once you have completed these steps, you should be able to open eclipse and find the use the files.
Msg me on Discord if you have any questions :smile::+1:
