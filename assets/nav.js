/**
 * Injects the shared navigation bar into #site-nav and marks the active link.
 * Call after DOM is ready: <script src="assets/nav.js" defer></script>
 */
(function () {
  const pages = [
    { href: 'index.html',   label: 'Overview' },
    { href: 'methods.html', label: 'Methods' },
    { href: 'results.html', label: 'Results' },
    { href: 'about.html',   label: 'About' },
  ];

  const current = window.location.pathname.split('/').pop() || 'index.html';

  const links = pages
    .map(({ href, label }) => {
      const isActive = href === current || (current === '' && href === 'index.html');
      return `<a href="${href}" class="${isActive ? 'active' : ''}">${label}</a>`;
    })
    .join('');

  const html = `
    <div style="max-width:1100px;margin:0 auto;padding:0 1.5rem;display:flex;align-items:center;justify-content:space-between;height:56px;position:relative;">
      <a href="index.html" class="nav-brand">STAD<span>-Seq</span></a>
      <button id="nav-toggle" aria-label="Toggle menu" aria-expanded="false">
        <span></span><span></span><span></span>
      </button>
      <nav id="nav-links" style="display:flex;align-items:center;gap:2rem;">
        ${links}
      </nav>
    </div>
  `;

  const container = document.getElementById('site-nav');
  if (container) {
    container.innerHTML = html;

    const toggle = document.getElementById('nav-toggle');
    const navLinks = document.getElementById('nav-links');

    toggle.addEventListener('click', () => {
      const open = navLinks.classList.toggle('open');
      toggle.setAttribute('aria-expanded', open);
    });

    // Close menu when a link is clicked (mobile)
    navLinks.querySelectorAll('a').forEach(a => {
      a.addEventListener('click', () => {
        navLinks.classList.remove('open');
        toggle.setAttribute('aria-expanded', 'false');
      });
    });
  }
})();
