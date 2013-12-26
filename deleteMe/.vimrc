"yao's vimrc
"Last modified: 20131115
set nocompatible
syntax on
set encoding=utf-8

" This allows re-use the same window when using multiple buffer, without
" saving first
set hidden

"Highlighten search result
"
set hlsearch

" Use case insensitive search, except when using capital letters
set ignorecase
set smartcase

" Display cursor position
set ruler

if has("autocmd")
  " Enable file type detection.
  " Use the default filetype settings, so that mail gets 'tw' set to 72,
  " 'cindent' is on in C files, etc.
  " Also load indent files, to automatically do language-dependent indenting
  filetype plugin indent on
endif

autocmd FileType python setlocal expandtab shiftwidth=4 softtabstop=4
autocmd FileType xml setlocal expandtab shiftwidth=2 softtabstop=2
autocmd FileType c,cpp,cc,h set number cindent expandtab shiftwidth=4 softtabstop=4

set guifont=Menlo\ Regular:h14
