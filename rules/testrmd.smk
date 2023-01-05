rule markdown:
    params:
        param1 = 'hello',
        param2 = 'world!'
    output:
        'markdown.html'
    shell:
        "Rscript -e \"rmarkdown::render('markdown.Rmd', params=list(param1='{params.param1}', param2='{params.param2}'))\""