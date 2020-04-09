#!/usr/bin/env nextflow

// testjes.nf



    print "hello"
    println "hello to you too"

    x="how are you"
    println x

    // strings
    a="hello there"
    print a + " ojo " + "\n"
    ojo="today"
    println '$ojo' // double quoted: support variables
    println "$ojo" // single quoted: do not!
    text = """   
    ojo ojoj \
    lalalal
    """ // multiline strings

    ''''
    myLongCmdline = """ blastp \
                    -in $input_query \
                    -out $output_file \
                    -db $blast_database \
                    -html
                    """

    result = myLongCmdline.execute().text
    '''


    //list
    lst=["good", "not so good", "I've been better"]
    println lst
    println lst[1]
    println lst.size()

    //map =? dictionary
    scores=["aswer21":5, "answer2":"good"]
    println scores
    println scores.answer2
    scores["answerX"]= 66
    println scores

    //conditional: if
    x=Math.random()
    if(x<0.5){
        println "loose"
    }
    else {
        println "woohoow"
        }






    // closure = can be argument to function
    square = { it * it }
    println square(9)

    array=[ 1, 2, 3, 4 ].collect(square)  //collect runs through each item in the array
    println array

    printMapClosure = { key, value ->
        println "$key = $value"
    }

    [ "Yue" : "Wu", "Mark" : "Williams", "Sudha" : "Kumari" ].each(printMapClosure)


myMap = ["China": 1 , "India" : 2, "USA" : 3]

result = 0
myMap.keySet().each( { result+= myMap[it] } )

println result


// regular expresion: ~/pattern/
''''
x1= assert 'foobar' =~/foo/  //return true
x2= assert 'foobar' ==~/foo/    // return FALSE
println x1 // does not work
println x2
'''

// string replacement
x = "colour".replaceFirst(/ou/, "o")
println x
// prints: color

y = "cheesecheese".replaceAll(/cheese/, "nice")
println y
// prints: nicenice

//capturing
programVersion = '2.7.3-beta'
(full, major, minor, patch, flavor) = (programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/)[0]

println full    // 2.7.3-beta
println major   // 2
println minor   // 7
println patch   // 3
println flavor  // beta

//removing part of string
// define the regexp pattern
wordStartsWithGr = ~/(?i)\s+Gr\w+/

// apply and verify the result
('Hello Groovy world!' - wordStartsWithGr) == 'Hello world!'
('Hi Grails users' - wordStartsWithGr) == 'Hi users'

println ('Remove first match of 5 letter word' - ~/\b\w{5}\b/) //remove  match of 5 letter word


// FILES
myFile = file('/home/hannelore/PROJECTHH/readme.txt')
//use glob pathµ
// use ** to also search beneath directory
listOfFiles = file('/home/hannelore/PROJECTHH/dataKP-assembly/KP-testset20/*.fna')
listOfFilesBIG = file('/home/hannelore/PROJECTHH/dataKP-assembly/*.fna')
''''
print myFile
myFile.text= "add this info"  //overwrite!
myFile.append= "readme" //add
myFile << 'Add a line more\n'
'''

// read line by line: dit werkt
file('/home/hannelore/PROJECTHH/Scripts/my_file.txt')
    .readLines()  //not for big files
    .each { println it }

//waarom werkt dit niet?
myFile2 = file('/home/hannelore/PROJECTHH/Scripts/my_file.txt')  //opm: does not create a file!
allelijn  = myFile2.readLines()
//!!! gebruik in en niet "":""
for( line in allelijn) {
    println line
}
//OPM: werkt niet voor append bvb!!???

//for big files: eachLine
count = 0
myFile2.eachLine {  str ->
        println "line ${count++}: $str"
    }

// advanced file reading/writing: X
''''

//create directories
myDir = file('anypath')
result = myDir.mkdir()
println result ? "OK" : "Cannot create directory: $myDir"

//create links: X

//copy and move files
myFile2.copyTo('new_name.txt')  //also with dirs, if not exist: created

myFile2.moveTo('anypath/new_file.txt') //also for directories and their content

//rename and delete
myFile3 = file('new_name.txt')
myFile3.renameTo('new_file_name.txt')

result = myFile3.delete()
println result ? "OK" : "Can delete: $myFile3"'

'''

//counting
def sample = file('new_name.txt')
println sample.countLines()

//voor multifastas
def sample = file('test.fasta')
println sample.countFasta()

// voor fastq
def sample = file('/data/sample.fastq')
println sample.countFastq()