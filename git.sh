rm -rf .git
git init
git add .
git commit -m "changes"
git branch -M main
git remote add origin https://github.com/PrantikDas/Cavity-Flow-with-Point-Particles
git push -u origin main -f
