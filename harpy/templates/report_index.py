from harpy import __version__

INDEXMD = """---
title: Welcome to My Landing Page
site:
  hide_outline: true
  hide_title_block: true
edit_url: null
---

{{button}}`Harpy Version {0} <https://github.com/pdimens/harpy>`

```{{image}} https://github.com/pdimens/harpy/blob/docs/static/logo_trans.png?raw=true
:class: col-page-right
:alt: The Harpy software logo
:width: 100%
:align: center
```

::::{{grid}} 1 1 2 3
:class: col-page-right
:::{{card}}
:header: Harpy Reports 📝
This is an aggregation of `.ipynb` reports produced by Harpy ([](https://doi.org/10.1093/bioadv/vbaf133)), rendered
in HTML by [MyST](https://mystmd.org/). Use the left sidebar
to navigate the directories and their reports.
:::

:::{{card}}
:header: Let us know of issues 🚩
If there are issues/errors in these reports, please [submit
an Issue](https://github.com/pdimens/harpy/issues/new/choose) on GitHub. 
:::

:::{{card}}
:header: What you will find 🔎
Stand-alone HTML reports created by other software (_e.g._ `fastp` or `MultiQC`) are not inlcluded into this aggregation due to
software limitations. 
:::
:::
""".format(__version__)