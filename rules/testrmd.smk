rule markdown:
    input:
        param1 = 'hello'
        param2 = 'world!'
    output:
        'markdown.html'
    shell:
        "Rscript -e \"rmarkdown::render('markdown.Rmd', params=list(param1='{input.param1}', param2='{input.param2}'))\""