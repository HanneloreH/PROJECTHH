



//Remove final "/"" from workingfolders and define them:
//function to trim final /
def trimFolder = { 
    it.endsWith("/") ? it[0..-2] : it
}
//do it for working folders


trainingx = trimFolder("$params.training")

println "folder is $trainingx"

