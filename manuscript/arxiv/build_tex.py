"""
Convert paper.md + supplementary_information.md -> arxiv/main.tex.

Pipeline:
  1. Pre-process markdown to convert author-year citations to \\cite{key}
     using a hand-curated map keyed off references.bib.
  2. Convert unicode glyphs to LaTeX.
  3. Replace markdown figure callouts with \\includegraphics blocks.
  4. Pipe through pandoc with a minimal LaTeX writer (no template).
  5. Wrap in a custom preamble (article class, natbib, graphicx, hyperref,
     amsmath, etc.).
  6. Append SI as appendix.

Run:  python build_tex.py
Output:  main.tex (in same directory)
"""

from __future__ import annotations
from pathlib import Path
import re
import subprocess

ARXIV = Path(__file__).resolve().parent
ROOT  = ARXIV.parent.parent           # repo root
MS    = ROOT / "manuscript"

# -------------------------------------------------------------------------
# 1. CITATION MAP — every inline author-year string we see in the body, and
#    the bib key(s) it resolves to. Keep exhaustive.
# -------------------------------------------------------------------------
CITE_MAP = {
    # Single-author / pair / et al.
    "Ang & Chen 2002":           "angchen2002",
    "Longin & Solnik 2001":      "longinsolnik2001",
    "Patton 2006":               "patton2006",
    "Daníelsson 2002":           "danielsson2002",
    "Daníelsson 2013":           "danielsson2013",
    "Brunnermeier & Pedersen 2009": "brunnermeier2009",
    "Ashby 1956":                "ashby1956",
    "Muhlert 2026":              "muhlert2026",
    "Perrow 1984":               "perrow1984",
    "Svolik 2019":               "svolik2019",
    "Levitsky & Ziblatt 2018":   "levitskyziblatt2018",
    "Chen et al. 2024":          "chen2024",
    "Sifat & Mohamad 2019":      "sifatmohamad2019",
    "Brogaard et al. 2018":      "brogaard2018",
    "McKay & Reis 2016":         "mckayreis2016",
    "IPBES 2019":                "ipbes2019",
    "Barnosky et al. 2011":      "barnosky2011",
    "Hänggi et al. 1990":        "hanggi1990",
    "Bacry et al. 2015":         "bacry2015",
    "Li et al. 2023":            "li2023alpacaeval",
    "Zou et al. 2023":           "zou2023universal",
    "Wei et al. 2023":           "wei2023jailbroken",
    "Mazeika et al. 2024":       "mazeika2024harmbench",
    "Sornette 2009":             "sornette2009",
    "Sornette & Ouillon 2012":   "sornette2012",
    "Bouchaud 2011":             "bouchaud2011",
    "Bouchaud & Potters 2003":   "bouchaud2003",
}

# Author-year continuation patterns — appear in text like
#   "Daníelsson [2002; 2013]"
# meaning two cites of the same author. Map these to compound bibkeys.
CONTINUATION_PATTERNS = [
    # (regex pattern that matches the bracketed continuation citation,
    #  list of bib keys to resolve to, given the immediately-preceding context)
]


# -------------------------------------------------------------------------
# 2. UNICODE -> LATEX glyph map
# -------------------------------------------------------------------------
GLYPHS = {
    # Greek
    "ρ": r"$\rho$",
    "μ": r"$\mu$",
    "σ": r"$\sigma$",
    "ξ": r"$\xi$",
    "α": r"$\alpha$",
    "β": r"$\beta$",
    "γ": r"$\gamma$",
    "η": r"$\eta$",
    "δ": r"$\delta$",
    "ε": r"$\varepsilon$",
    "τ": r"$\tau$",
    "Φ": r"$\Phi$",
    "Δ": r"$\Delta$",
    # Math operators
    "²": r"$^2$",
    "³": r"$^3$",
    "⁻¹": r"$^{-1}$",
    "⁻²": r"$^{-2}$",
    "⁻³": r"$^{-3}$",
    "⁻⁴": r"$^{-4}$",
    "⁻⁵": r"$^{-5}$",
    "⁻⁶": r"$^{-6}$",
    "⁻⁷": r"$^{-7}$",
    "⁻¹⁵": r"$^{-15}$",
    "ⁿ": r"$^n$",
    "⁰": r"$^0$",
    "¹": r"$^1$",
    "⁴": r"$^4$",
    "⁵": r"$^5$",
    "⁶": r"$^6$",
    "⁷": r"$^7$",
    "⁸": r"$^8$",
    "⁹": r"$^9$",
    "ₑ": r"$_\text{e}$",
    "ᶜ": r"$^c$",
    "√": r"$\sqrt{}$",
    "×": r"$\times$",
    "−": "-",
    "→": r"$\to$",
    "≤": r"$\leq$",
    "≥": r"$\geq$",
    "≈": r"$\approx$",
    "∈": r"$\in$",
    "∞": r"$\infty$",
    "·": r"$\cdot$",
    # Quotes / dashes (rendered cleanly by LaTeX itself if we use right form)
    "–": "--",
    "—": "---",
    "“": "``",   # left double quote
    "”": "''",   # right double quote
    "‘": "`",    # left single quote
    "’": "'",    # right single quote
    "…": r"\ldots ",
    "°": r"$^\circ$",
    "ℓ": r"$\ell$",
}

# Special LaTeX characters needing escape OUTSIDE math mode.
SPECIAL_LATEX = {
    "%": r"\%",
    "&": r"\&",
    "#": r"\#",
    "_": r"\_",
}


def escape_latex_outside_math(text: str) -> str:
    """Escape %, &, #, _ everywhere except inside $...$ math. Pandoc handles
    most of this for us in body text; we only need to fix raw-pasted spans."""
    return text  # pandoc handles this; no-op


def replace_glyphs(text: str) -> str:
    # Order matters: longer multi-char unicode FIRST (e.g. ⁻¹⁵ before ⁻¹)
    for glyph, latex in sorted(GLYPHS.items(), key=lambda kv: -len(kv[0])):
        text = text.replace(glyph, latex)
    return text


def replace_citations(text: str) -> str:
    """Replace bracketed author-year citations with \\cite{} commands.

    Strategy:
      - For every bracketed expression containing a 4-digit year, parse the
        contents and produce the corresponding \\cite{key1,key2,...}.
      - Continuation citations like "[2002; 2013]" are converted by walking
        backwards to the most recent author name and resolving via the map.
    """

    # Step A: Inline bracketed cites with full author-year inside the brackets.
    def repl(match):
        body = match.group(1).strip()
        # Strip embedded newlines that pandoc-parsed markdown sometimes inserts
        body = re.sub(r"\s+", " ", body)
        # Split on ;
        parts = [p.strip() for p in body.split(";")]
        keys = []
        for p in parts:
            # Try direct lookup
            if p in CITE_MAP:
                keys.append(CITE_MAP[p])
                continue
            # Try with "et al."
            # Try matching author + year pattern; if it matches a known map entry, use it
            found = False
            for k_pattern, bibkey in CITE_MAP.items():
                if k_pattern == p:
                    keys.append(bibkey)
                    found = True
                    break
            if not found:
                # Year-only continuation -- mark for second pass
                m_year = re.match(r"^(\d{4})$", p)
                if m_year:
                    # Defer: leave as a sentinel
                    keys.append(f"__YEAR_{p}__")
                else:
                    # Unknown citation pattern -- leave verbatim with warning marker
                    keys.append(f"__UNKNOWN_{p.replace(' ', '_')}__")
        return r"\cite{" + ",".join(keys) + r"}"

    text = re.sub(r"\[((?:[^\[\]\n]|\n[ \t]*[^\[\]\n])*?\b\d{4}\b(?:[^\[\]]|\n[ \t]*[^\[\]\n])*?)\]", repl, text)

    # Step B: Resolve __YEAR_2002__ continuations by walking backwards in text
    #          to find the most recent author name. The corresponding
    #          \cite{} replacement produces a key like "danielsson2002".
    AUTHOR_TO_PREFIX = {
        "Daníelsson": "danielsson",
        "Sornette":   "sornette",
        "Bouchaud":   "bouchaud",
        "Brunnermeier": "brunnermeier",
        "Ashby":      "ashby",
        "Perrow":     "perrow",
        "Svolik":     "svolik",
        "Patton":     "patton",
        "Hänggi":     "hanggi",
        "IPBES":      "ipbes",
        "Zou":        "zou2023universal_",   # year is appended; for Zou map to fixed key
        "Wei":        "wei2023jailbroken_",
        "Mazeika":    "mazeika2024harmbench_",
        "Li":         "li2023alpacaeval_",
        "Bacry":      "bacry",
    }
    # Some bib keys do not follow the simple "<prefix><year>" convention.
    # Override the resolved key here:
    KEY_OVERRIDES = {
        "zou2023universal_2023":     "zou2023universal",
        "wei2023jailbroken_2023":    "wei2023jailbroken",
        "mazeika2024harmbench_2024": "mazeika2024harmbench",
        "li2023alpacaeval_2023":     "li2023alpacaeval",
    }

    def resolve_year_sentinels(text: str) -> str:
        out = []
        i = 0
        while i < len(text):
            m = re.search(r"__YEAR_(\d{4})__", text[i:])
            if not m:
                out.append(text[i:])
                break
            start = i + m.start()
            year = m.group(1)
            # Walk backwards from `start` to find the most recent author name
            preceding = text[:start]
            best_author = None
            best_pos = -1
            for author, prefix in AUTHOR_TO_PREFIX.items():
                pos = preceding.rfind(author)
                if pos > best_pos:
                    best_pos = pos
                    best_author = (author, prefix)
            if best_author is None:
                replacement = "??UNKNOWNAUTHOR" + year
            else:
                replacement = best_author[1] + year
                # Apply override if the synthesized key isn't actually in the bib
                if replacement in KEY_OVERRIDES:
                    replacement = KEY_OVERRIDES[replacement]
            out.append(text[i:start])
            out.append(replacement)
            i = start + m.end() - m.start()
        return "".join(out)

    text = resolve_year_sentinels(text)

    return text


def wrap_figures_in_latex(tex: str, label_prefix: str = "fig",
                          path_prefix: str = "figures/",
                          si_prefix: bool = False) -> str:
    """After pandoc has converted markdown to LaTeX, recognise the converted
    figure-caption pattern  \\textbf{Figure N.} \\emph{title.} body text...
    and wrap into a proper LaTeX figure environment with \\caption and
    \\label.

    The body text continues until the next \\textbf{Figure ...} marker or
    end of document.
    """
    # Stop the caption at: next figure marker; OR a section/subsection heading;
    # OR a pandoc horizontal rule; OR end of input.
    STOP_LOOKAHEAD = (
        r"(?=\\textbf\{Figure S?\d+\.\}"
        r"|\\(?:section|subsection|subsubsection)\*?\{"
        r"|\\begin\{center\}\\rule"
        r"|\Z)"
    )
    if si_prefix:
        pattern = re.compile(
            r"\\textbf\{Figure S(\d+)\.\}\s*\\emph\{([^}]+?)\}\s*(.*?)" + STOP_LOOKAHEAD,
            re.DOTALL,
        )
    else:
        pattern = re.compile(
            r"\\textbf\{Figure (\d+)\.\}\s*\\emph\{([^}]+?)\}\s*(.*?)" + STOP_LOOKAHEAD,
            re.DOTALL,
        )

    def repl(m):
        n = m.group(1)
        title = m.group(2)
        body = m.group(3).strip()
        # Collapse internal blank lines (caption can't contain paragraph breaks)
        body = re.sub(r"\n\s*\n+", " ", body)
        body = re.sub(r"\s+", " ", body).strip()
        title = title.strip()
        # Build the figure environment
        prefix = "figS" if si_prefix else "fig"
        return (
            "\n\n\\begin{figure*}[!htbp]\n"
            "\\centering\n"
            f"\\includegraphics[width=0.95\\textwidth]{{{path_prefix}{prefix}{n}.png}}\n"
            f"\\caption{{\\textit{{{title}}} {body}}}\n"
            f"\\label{{fig:{prefix}{n}}}\n"
            "\\end{figure*}\n\n"
        )

    return pattern.sub(repl, tex)


def preprocess(md_text: str, is_si: bool = False) -> str:
    text = md_text

    # 0. Strip the leading H1 title — for the main paper it's already in
    #    \title{} via \maketitle; for the SI it's in our appendix wrapper.
    text = re.sub(r"^# [^\n]*\n", "", text, count=1)

    if is_si:
        # SI section/subsection headings include "§S1" (sections) and "S1.1"
        # (subsections). The appendix counter would re-prefix them. Strip both
        # forms.
        text = re.sub(r"(^#+\s+)§?S\d+(?:\.\d+)?\s+", r"\1", text, flags=re.MULTILINE)
        # Strip the "References for the Supplementary Information" section —
        # natbib's \bibliography already prints all SI-cited refs, so the
        # manual list is redundant and would otherwise be numbered as S12,
        # contradicting the cover letter's "§S1–§S11".
        text = re.sub(r"^## References for the Supplementary Information.*?(?=\n## |\Z)",
                      "", text, count=1, flags=re.MULTILINE | re.DOTALL)

    if not is_si:
        # Strip the "## Figures" heading from the main paper — the figure
        # environments are already inserted inline via wrap_figures_in_latex,
        # so this heading would render as an empty section.
        text = re.sub(r"^## Figures\s*\n", "", text, flags=re.MULTILINE)

    # 1. Citations
    text = replace_citations(text)

    # 2. Figure environments are handled POST-pandoc (so pandoc converts
    #    markdown formatting inside captions like *italic*, **bold**, ρ, μ etc).

    # 3. Unicode glyphs are handled in the LaTeX preamble via newunicodechar.

    return text


def run_pandoc(md_text: str) -> str:
    """Convert preprocessed markdown to LaTeX body via pandoc."""
    proc = subprocess.run(
        ["pandoc",
         "--from=markdown+raw_tex",
         "--to=latex",
         "--wrap=none",
         "--no-highlight",
         # Demote everything by one level: ## becomes \section, ### becomes \subsection.
         # The H1 (# Title) is already supplied via \title{} + \maketitle;
         # we strip it from the source before sending to pandoc.
         "--shift-heading-level-by=-1"],
        input=md_text.encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return proc.stdout.decode("utf-8")


def postprocess_latex(tex: str) -> str:
    """After pandoc emits the LaTeX, fix up:

      1. Display math \\[ ... \\tag{N} ... \\] -> equation environment
         (pandoc converts $$...$$ to \\[ ... \\] but \\tag{} requires
         equation/equation*).
    """
    # \[ ... \tag{...} ... \] -> \begin{equation*} ... \end{equation*}
    # IMPORTANT: collapse the equation body to a single line — blank lines
    # inside math trigger paragraph breaks which are illegal in math mode.
    pattern = re.compile(r"\\\[(.*?)\\\]", re.DOTALL)

    def repl(m):
        body = m.group(1)
        if r"\tag{" in body:
            # collapse all whitespace runs (including newlines) to single spaces
            collapsed = re.sub(r"\s+", " ", body).strip()
            return "\n\\begin{equation*}\n" + collapsed + "\n\\end{equation*}\n"
        # Even for non-tag display math, ensure no blank lines inside \[ \]
        collapsed = re.sub(r"\n\s*\n+", "\n", body)
        return "\\[\n" + collapsed.strip() + "\n\\]"

    tex = pattern.sub(repl, tex)

    # Convert \texttt{path/with/slashes_or_dots} to \protect\path|...|, which
    # is the url-package macro for paths: typewriter font, breakable at
    # / _ . - chars. \protect makes it safe inside \caption{} (moving args).
    # Only target patterns containing at least one / (a clear file-path
    # signal); short bits like \texttt{*.csv} or \texttt{s_0} stay as-is.
    def path_repl(m):
        body = m.group(1)
        # \path uses verbatim semantics — undo pandoc's backslash-escapes.
        body = body.replace(r"\_", "_").replace(r"\&", "&").replace(r"\%", "%")
        # Choose a delimiter character not present in the path body.
        for delim in "|!@^+~":
            if delim not in body:
                return f"\\protect\\path{delim}{body}{delim}"
        return m.group(0)

    # Path-shaped: contains a / AND only ASCII path-safe chars
    # (letters, digits, /, _, ., -, *, the LaTeX-escape \\, &, %).
    PATH_BODY = r"[A-Za-z0-9_./\-*\\&%]*?/[A-Za-z0-9_./\-*\\&%]*?"
    tex = re.sub(r"\\texttt\{(" + PATH_BODY + r")\}", path_repl, tex)
    return tex


def build_main_tex():
    paper_md = (MS / "paper.md").read_text(encoding="utf-8")
    si_md    = (MS / "supplementary_information.md").read_text(encoding="utf-8")

    paper_md_pre = preprocess(paper_md, is_si=False)
    si_md_pre    = preprocess(si_md, is_si=True)

    paper_tex = wrap_figures_in_latex(postprocess_latex(run_pandoc(paper_md_pre)),
                                       si_prefix=False)
    si_tex    = wrap_figures_in_latex(postprocess_latex(run_pandoc(si_md_pre)),
                                       si_prefix=True)

    # Strip pandoc's outer document boilerplate if any
    # (with no template, pandoc gives us only the body — fine.)

    preamble = r"""\documentclass[11pt,letterpaper]{article}

% Compile with: xelatex main.tex; bibtex main; xelatex main.tex; xelatex main.tex
% xelatex natively understands UTF-8, no inputenc needed.

\usepackage[margin=1in]{geometry}
\usepackage{fontspec}            % xelatex font selection
\usepackage{amsmath,amssymb,amsthm}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage[numbers,sort&compress]{natbib}
\usepackage{booktabs}
\usepackage{caption}
\usepackage{longtable}
\usepackage{array}
\usepackage{calc}
% pandoc emits \real{0.2} etc. inside table column specs; provide it.
\providecommand{\real}[1]{#1}
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
% pandoc emits \def\LTcaptype{none} to suppress table captions; define the
% counter so longtable doesn't error.
\newcounter{none}
\usepackage{microtype}
\usepackage{titlesec}
\usepackage{newunicodechar}
\usepackage{url}
% Allow line-breaks inside \url{} / \path{} at slashes, underscores,
% dots, hyphens. Output keeps typewriter font.
\PassOptionsToPackage{hyphens}{url}
% No stretchy spacing at break points — keeps path characters tight when
% the line is justified. Breaks are still enabled by \UrlBreaks below.
\Urlmuskip=0mu\relax
% Add break-points after / _ . - so file paths split cleanly.
\renewcommand{\UrlBreaks}{\do\/\do\_\do\.\do\-\do\&\do\?\do\=\do\,\do\;}
% Loosen TeX's line-breaking globally, to absorb any remaining overflow.
\sloppy
\setlength{\emergencystretch}{3em}

% Unicode -> math/text mappings (so utf8 source compiles with pdflatex)
% Greek letters
\newunicodechar{ρ}{\ensuremath{\rho}}
\newunicodechar{μ}{\ensuremath{\mu}}
\newunicodechar{σ}{\ensuremath{\sigma}}
\newunicodechar{ξ}{\ensuremath{\xi}}
\newunicodechar{α}{\ensuremath{\alpha}}
\newunicodechar{β}{\ensuremath{\beta}}
\newunicodechar{γ}{\ensuremath{\gamma}}
\newunicodechar{η}{\ensuremath{\eta}}
\newunicodechar{δ}{\ensuremath{\delta}}
\newunicodechar{ε}{\ensuremath{\varepsilon}}
\newunicodechar{τ}{\ensuremath{\tau}}
\newunicodechar{π}{\ensuremath{\pi}}
\newunicodechar{Φ}{\ensuremath{\Phi}}
\newunicodechar{Δ}{\ensuremath{\Delta}}
% Math operators / scripts
\newunicodechar{²}{\ensuremath{^{2}}}
\newunicodechar{³}{\ensuremath{^{3}}}
\newunicodechar{⁰}{\ensuremath{^{0}}}
\newunicodechar{¹}{\ensuremath{^{1}}}
\newunicodechar{⁴}{\ensuremath{^{4}}}
\newunicodechar{⁵}{\ensuremath{^{5}}}
\newunicodechar{⁶}{\ensuremath{^{6}}}
\newunicodechar{⁷}{\ensuremath{^{7}}}
\newunicodechar{⁸}{\ensuremath{^{8}}}
\newunicodechar{⁹}{\ensuremath{^{9}}}
\newunicodechar{ⁿ}{\ensuremath{^{n}}}
\newunicodechar{⁻}{\ensuremath{^{-}}}
\newunicodechar{ₑ}{\ensuremath{_{\mathrm{e}}}}
\newunicodechar{ᶜ}{\ensuremath{^{c}}}
\newunicodechar{₀}{\ensuremath{_{0}}}
\newunicodechar{₊}{\ensuremath{_{+}}}
\newunicodechar{₋}{\ensuremath{_{-}}}
\newunicodechar{√}{\ensuremath{\sqrt{\,}}}
\newunicodechar{×}{\ensuremath{\times}}
\newunicodechar{−}{\ensuremath{-}}
\newunicodechar{→}{\ensuremath{\rightarrow}}
\newunicodechar{≤}{\ensuremath{\leq}}
\newunicodechar{≥}{\ensuremath{\geq}}
\newunicodechar{≈}{\ensuremath{\approx}}
\newunicodechar{∈}{\ensuremath{\in}}
\newunicodechar{∞}{\ensuremath{\infty}}
\newunicodechar{·}{\ensuremath{\cdot}}
\newunicodechar{±}{\ensuremath{\pm}}
\newunicodechar{∝}{\ensuremath{\propto}}
\newunicodechar{∧}{\ensuremath{\wedge}}
\newunicodechar{∨}{\ensuremath{\vee}}
\newunicodechar{∑}{\ensuremath{\sum}}
\newunicodechar{≠}{\ensuremath{\neq}}
\newunicodechar{≫}{\ensuremath{\gg}}
\newunicodechar{≪}{\ensuremath{\ll}}
% Subscript Latin letters (used in symbols like h_max, J_eff, c_log)
\newunicodechar{ₐ}{\ensuremath{_{a}}}
\newunicodechar{ₓ}{\ensuremath{_{x}}}
\newunicodechar{ₖ}{\ensuremath{_{k}}}
\newunicodechar{ₘ}{\ensuremath{_{m}}}
\newunicodechar{ₚ}{\ensuremath{_{p}}}
\newunicodechar{ₛ}{\ensuremath{_{s}}}
\newunicodechar{ₜ}{\ensuremath{_{t}}}
\newunicodechar{ⱼ}{\ensuremath{_{j}}}
\newunicodechar{≳}{\ensuremath{\gtrsim}}
\newunicodechar{≲}{\ensuremath{\lesssim}}
\newunicodechar{⌊}{\ensuremath{\lfloor}}
\newunicodechar{⌋}{\ensuremath{\rfloor}}
\newunicodechar{⟨}{\ensuremath{\langle}}
\newunicodechar{⟩}{\ensuremath{\rangle}}
\newunicodechar{ℓ}{\ensuremath{\ell}}
% Punctuation
\newunicodechar{§}{\S{}}
\newunicodechar{°}{\textdegree}
\newunicodechar{…}{\ldots}
% Phonetic / IPA / specialty (best-effort fallbacks)
\newunicodechar{ˡ}{\ensuremath{^{l}}}
\newunicodechar{ˣ}{\ensuremath{^{x}}}
\newunicodechar{̄}{}
\newunicodechar{ᐟ}{}
\newunicodechar{ᴅ}{D}
\newunicodechar{ᵒ}{\ensuremath{^{o}}}
\newunicodechar{ᵢ}{\ensuremath{_{i}}}
\newunicodechar{ᵦ}{\ensuremath{_{\beta}}}
\newunicodechar{ő}{\H{o}}

\hypersetup{
  colorlinks=true,
  linkcolor=blue!50!black,
  citecolor=blue!50!black,
  urlcolor=blue!50!black,
}

\titleformat*{\section}{\Large\bfseries}
\titleformat*{\subsection}{\large\bfseries}
\titleformat*{\subsubsection}{\normalsize\bfseries\itshape}

\setlength{\parskip}{4pt}
\setlength{\parindent}{0pt}

\title{Why Fixed Protections Fail Under Rising Coordination:\\
A Structural Failure-Mode Transition in Coupled Systems}
\author{Chowon Jung\\
\small\texttt{dancing4am@gmail.com}}
\date{\today}

\begin{document}
\maketitle

"""

    middle = r"""

\bibliographystyle{unsrtnat}
\bibliography{references}

\appendix
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thesubsection}{S\arabic{section}.\arabic{subsection}}
\renewcommand{\thefigure}{S\arabic{figure}}
\setcounter{figure}{0}
\setcounter{section}{0}

\section*{Supplementary Information}

"""

    end = r"""

\end{document}
"""

    # Combined manuscript (paper + SI as appendix) — for arXiv.
    out_combined = preamble + paper_tex + middle + si_tex + end

    # Paper-only — for NatComms "main manuscript" upload field.
    paper_only_tail = r"""

\bibliographystyle{unsrtnat}
\bibliography{references}

\end{document}
"""
    out_paper_only = preamble + paper_tex + paper_only_tail

    # SI-only — for NatComms "supplementary information" upload field.
    si_preamble = preamble.replace(
        r"\title{Why Fixed Protections Fail Under Rising Coordination:\\"
        "\nA Structural Failure-Mode Transition in Coupled Systems}",
        r"\title{Supplementary Information \\ \large Why Fixed Protections "
        r"Fail Under Rising Coordination:\\ A Structural Failure-Mode "
        r"Transition in Coupled Systems}"
    )
    out_si_only = (
        si_preamble
        + r"""
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thesubsection}{S\arabic{section}.\arabic{subsection}}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\theequation}{S\arabic{equation}}
\setcounter{figure}{0}
\setcounter{section}{0}
\setcounter{equation}{0}

"""
        + si_tex
        + r"""

\bibliographystyle{unsrtnat}
\bibliography{references}

\end{document}
"""
    )

    # Normalize line endings: Pandoc on Windows emits \r\n; Python's text-mode
    # write would then add another \r, yielding \r\r\n which TeX reads as a
    # paragraph break. Force LF-only with newline="" + binary write.
    for fn, content in [
        ("main.tex",          out_combined),
        ("paper_only.tex",    out_paper_only),
        ("supplementary.tex", out_si_only),
    ]:
        normalized = content.replace("\r\n", "\n").replace("\r", "\n")
        (ARXIV / fn).write_bytes(normalized.encode("utf-8"))
        print(f"wrote {ARXIV / fn}: {len(normalized):,} chars (LF-only)")
    out = out_combined  # for sanity-check section below

    # Sanity-check: any unresolved citations or glyphs?
    unresolved_cites = re.findall(r"__UNKNOWN_[^\s\\]+", out)
    unresolved_years = re.findall(r"__YEAR_\d+__", out)
    unicode_left = []
    for ch in out:
        if ord(ch) > 127 and ch not in "‘’“”−–—…":
            unicode_left.append(ch)
    print(f"  unresolved citation patterns: {len(set(unresolved_cites))}")
    for u in sorted(set(unresolved_cites))[:10]:
        print(f"    {u}")
    print(f"  unresolved year continuations: {len(set(unresolved_years))}")
    print(f"  remaining non-ASCII chars: {len(set(unicode_left))}")
    for u in sorted(set(unicode_left))[:20]:
        print(f"    {repr(u)} = U+{ord(u):04X}")


if __name__ == "__main__":
    build_main_tex()
