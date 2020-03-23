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

    """
    #!/usr/bin/perl

    print 'Hi there!' . '\n';
    """

}

process pyStuff {

    """
    #!/usr/bin/python

    x = 'Hello'
    y = 'world!'
    print "%s - %s" % (x,y)
    """

}