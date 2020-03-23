#!/usr/bin/env nextflow

process sayHello {
    """
    echo 'Hello world!' > file
    """
}

///opm: de m """ vooral nodig bij multi-strings!

'''
process 5blocks {

   [ directives ]

   input:
    < process inputs >

   output:
    < process outputs >

   when:
    < condition >

   [script|shell|exec]:
   < user script to be executed >

}
'''

PATH="test"
// use single quotes! om variable van BASH (en niet van de pipeline) te gebruiken
process printPath2 {
  output:
  stdout result
   """
 echo The path is: $PATH
   """
}

process printPath {
   output:
  stdout result2
   '''
   echo The path is: $PATH
  '''
}


//andere print: nu printf

// om resultaten te zien
result.view { it }
result2.view { it }


// other languages
process perlStuff {
  output:
  stdout result3
    """
    #!/usr/bin/env perl

    print 'Hi there!' . '\n';
    """
}
/// tip: Since the actual location of the interpreter binary file can change across platforms, to make your scripts more portable it is wise to use the env shell command followed by the interpreter's name, instead of the absolute path of it. Thus, the shebang declaration for a Perl script, for example, would look like: #!/usr/bin/env perl, does not work for python(?)

process pyStuff {
  output:
  stdout result4
    """
    #!/usr/bin/python

    x = 'Hello'
    y = 'world!'
    print "%s - %s" % (x,y)
    """
}

result3.view { it }
result4.view { it }

// conditional scripts

'''
process align {
    input:
    file seq_to_aln from sequences

    script:
    if( mode == 'tcoffee' )
        """
        t_coffee -in $seq_to_aln > out_file
        """

    else if( mode == 'mafft' )
        """
        mafft --anysymbol --parttree --quiet $seq_to_aln > out_file
        """

    else if( mode == 'clustalo' )
        """
        clustalo -i $seq_to_aln -o out_file
        """

    else
        error "Invalid alignment mode: ${mode}"

}
'''

//templates
// opm: my_script.sh moet in folder templates zitten in same dir of path geven
process template_example {

    input:
    val STR from 'this', 'that'

    output:
    stdout result5

    script:
    template 'my_script.sh'
}
result5.view { it }

//shell
process myTask {
    input:
    val str from 'Hello', 'Hola', 'Bonjour'

    output:
    stdout result6

    shell:
    '''
    echo User $USER says !{str}
    '''
}
result6.view { it }

//vs script
process myTask2 {
    input:
    val str2 from 'Hello', 'Hola', 'Bonjour'

    output:
    stdout result8

    script:
    """
    echo User \$USER says $str2
    """
}
result8.view { it }


// vs exec
x = Channel.from( 'a', 'b', 'c')
process simpleSum {
    input:
    val x
    output:
    stdout result9
    exec:
    println "Hello Mr. $x"
}
result9.view { it } // werkt niet???

//input:
  //<input qualifier> <input name> [from <source channel>] [attributes]
  //qualifiers: val (uit script), env (set var), file, path, stdin, tuple, each

  num = Channel.from( 1, 2, 3 )

process basicExample {
  input:
  val x from num // = val num
  output:
  stdout result10
  "echo process job $x"

}
result10.view { it }

'''
// maar waar is result (nog geen output)
proteins = Channel.fromPath( 'fnas/*.fna' )

process blastThemAll {
  input:
  file query_file from proteins //=file proteins

  "blastp -query ${query_file} -db nr"  //=$proteins
}'''

