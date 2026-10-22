// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import starlightUiTweaks from 'starlight-ui-tweaks'
import starlightScrollToTop from 'starlight-scroll-to-top'
import codeblocks from 'starlight-codeblocks';
import mermaid from 'astro-mermaid';

export default defineConfig({
  site: 'https://blinkseq.github.io/harpy',
  base: '/harpy',
	integrations: [
    mermaid(),
    starlight({
      customCss: ['./src/style/global.css'],
      title: '4.3',
      logo: {
        light: '/src/assets/logo_trans.png',
        dark: '/src/assets/logo_trans.png',
      },
      favicon: '/favicon.png',
			social: [{ icon: 'github', label: 'GitHub', href: 'https://github.com/BLinkseq/harpy' }],
			lastUpdated: true,
      pagination: false,
      editLink: {
        baseUrl: 'https://github.com/blinkseq/harpy/edit/main/docs/',
      },
      plugins: [
        codeblocks(),
        starlightUiTweaks({
        navbarLinks: [
           { label: "BLink-seq", href: "https://blinkseq.github.io"},
        ],
      }),
        starlightScrollToTop({
          showProgressRing: true,
          progressRingColor: '#83b4eb',
          borderRadius: '50',
        }),
      ],
		}),
	],
});
