




<!DOCTYPE html>
<html>
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# object: http://ogp.me/ns/object# article: http://ogp.me/ns/article# profile: http://ogp.me/ns/profile#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <title>GUI/Source/Processors/Editors/ElectrodeButtons.cpp at master · ankit--sethi/GUI</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <meta property="fb:app_id" content="1401488693436528"/>

      <meta content="@github" name="twitter:site" /><meta content="summary" name="twitter:card" /><meta content="ankit--sethi/GUI" name="twitter:title" /><meta content="GUI - Software for electrophysiology data acquisition" name="twitter:description" /><meta content="https://2.gravatar.com/avatar/b17a0b715c09bb2cd11e4abc02dd1663?d=https%3A%2F%2Fidenticons.github.com%2Fa0b87afc0689055d878e66270631162e.png&amp;r=x&amp;s=400" name="twitter:image:src" />
<meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="https://2.gravatar.com/avatar/b17a0b715c09bb2cd11e4abc02dd1663?d=https%3A%2F%2Fidenticons.github.com%2Fa0b87afc0689055d878e66270631162e.png&amp;r=x&amp;s=400" property="og:image" /><meta content="ankit--sethi/GUI" property="og:title" /><meta content="https://github.com/ankit--sethi/GUI" property="og:url" /><meta content="GUI - Software for electrophysiology data acquisition" property="og:description" />

    <meta name="hostname" content="github-fe122-cp1-prd.iad.github.net">
    <meta name="ruby" content="ruby 2.1.0p0-github-tcmalloc (87d8860372) [x86_64-linux]">
    <link rel="assets" href="https://github.global.ssl.fastly.net/">
    <link rel="conduit-xhr" href="https://ghconduit.com:25035/">
    <link rel="xhr-socket" href="/_sockets" />


    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="802AA11E:22DE:9C5BF4:52FAA786" name="octolytics-dimension-request_id" /><meta content="4641264" name="octolytics-actor-id" /><meta content="ankit--sethi" name="octolytics-actor-login" /><meta content="a21d1a93983dfd6b36a8b9a4d0ce1a96d2cb7e1c55635aba5fb551193a14bfc5" name="octolytics-actor-hash" />
    

    
    
    <link rel="icon" type="image/x-icon" href="/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="wDXFYi/4cJYkfprIO/l9R4DNrLucl0avlTBH86JJUgQ=" name="csrf-token" />

    <link href="https://github.global.ssl.fastly.net/assets/github-112d5f11fb1e8b126cb566609f8612a2b8a1dacc.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://github.global.ssl.fastly.net/assets/github2-da9d3c18bdc47da61291dd716c6cab8877cbdd70.css" media="all" rel="stylesheet" type="text/css" />
    


      <script src="https://github.global.ssl.fastly.net/assets/frameworks-752c70f2b89dcf2d1f948637afa35a3285fe6424.js" type="text/javascript"></script>
      <script async="async" defer="defer" src="https://github.global.ssl.fastly.net/assets/github-6ecfdea6c072948b639306b2bc266080946acf10.js" type="text/javascript"></script>
      
      <meta http-equiv="x-pjax-version" content="845eb3035b1f466b59548dbccbdaded3">

        <link data-pjax-transient rel='permalink' href='/ankit--sethi/GUI/blob/607c7ac003e80f8bf094a655a13a638c8247cf72/Source/Processors/Editors/ElectrodeButtons.cpp'>

  <meta name="description" content="GUI - Software for electrophysiology data acquisition" />

  <meta content="4641264" name="octolytics-dimension-user_id" /><meta content="ankit--sethi" name="octolytics-dimension-user_login" /><meta content="10552631" name="octolytics-dimension-repository_id" /><meta content="ankit--sethi/GUI" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="true" name="octolytics-dimension-repository_is_fork" /><meta content="3349728" name="octolytics-dimension-repository_parent_id" /><meta content="open-ephys/GUI" name="octolytics-dimension-repository_parent_nwo" /><meta content="3349728" name="octolytics-dimension-repository_network_root_id" /><meta content="open-ephys/GUI" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/ankit--sethi/GUI/commits/master.atom" rel="alternate" title="Recent Commits to GUI:master" type="application/atom+xml" />

  </head>


  <body class="logged_in  env-production linux vis-public fork page-blob">
    <div class="wrapper">
      
      
      
      


      <div class="header header-logged-in true">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/">
  <span class="mega-octicon octicon-mark-github"></span>
</a>

    
    <a href="/notifications" class="notification-indicator tooltipped downwards" data-gotokey="n" title="You have no unread notifications">
        <span class="mail-status all-read"></span>
</a>

      <div class="command-bar js-command-bar  in-repository">
          <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<input type="text" data-hotkey=" s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    data-username="ankit--sethi"
      data-repo="ankit--sethi/GUI"
      data-branch="master"
      data-sha="d8899fb0735ec99b7f9c8599971c648e2d310ad0"
  >

    <input type="hidden" name="nwo" value="ankit--sethi/GUI" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="octicon help tooltipped downwards" title="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
        <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
            <li><a href="https://gist.github.com">Gist</a></li>
            <li><a href="/blog">Blog</a></li>
          <li><a href="https://help.github.com">Help</a></li>
        </ul>
      </div>

    


  <ul id="user-links">
    <li>
      <a href="/ankit--sethi" class="name">
        <img alt="ankit--sethi" class=" js-avatar" data-user="4641264" height="20" src="https://1.gravatar.com/avatar/b17a0b715c09bb2cd11e4abc02dd1663?d=https%3A%2F%2Fidenticons.github.com%2Fa0b87afc0689055d878e66270631162e.png&amp;r=x&amp;s=140" width="20" /> ankit--sethi
      </a>
    </li>

    <li class="new-menu dropdown-toggle js-menu-container">
      <a href="#" class="js-menu-target tooltipped downwards" title="Create new..." aria-label="Create new...">
        <span class="octicon octicon-plus"></span>
        <span class="dropdown-arrow"></span>
      </a>

      <div class="js-menu-content">
      </div>
    </li>

    <li>
      <a href="/settings/profile" id="account_settings"
        class="tooltipped downwards"
        aria-label="Account settings "
        title="Account settings ">
        <span class="octicon octicon-tools"></span>
      </a>
    </li>
    <li>
      <a class="tooltipped downwards" href="/logout" data-method="post" id="logout" title="Sign out" aria-label="Sign out">
        <span class="octicon octicon-log-out"></span>
      </a>
    </li>

  </ul>

<div class="js-new-dropdown-contents hidden">
  

<ul class="dropdown-menu">
  <li>
    <a href="/new"><span class="octicon octicon-repo-create"></span> New repository</a>
  </li>
  <li>
    <a href="/organizations/new"><span class="octicon octicon-organization"></span> New organization</a>
  </li>


    <li class="section-title">
      <span title="ankit--sethi/GUI">This repository</span>
    </li>
      <li>
        <a href="/ankit--sethi/GUI/settings/collaboration"><span class="octicon octicon-person-add"></span> New collaborator</a>
      </li>
</ul>

</div>


    
  </div>
</div>

      

      




          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">

    <li class="subscription">
      <form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="wDXFYi/4cJYkfprIO/l9R4DNrLucl0avlTBH86JJUgQ=" /></div>  <input id="repository_id" name="repository_id" type="hidden" value="10552631" />

    <div class="select-menu js-menu-container js-select-menu">
      <a class="social-count js-social-count" href="/ankit--sethi/GUI/watchers">
        1
      </a>
      <span class="minibutton select-menu-button with-count js-menu-target" role="button" tabindex="0">
        <span class="js-select-button">
          <span class="octicon octicon-eye-unwatch"></span>
          Unwatch
        </span>
      </span>

      <div class="select-menu-modal-holder">
        <div class="select-menu-modal subscription-menu-modal js-menu-content">
          <div class="select-menu-header">
            <span class="select-menu-title">Notification status</span>
            <span class="octicon octicon-remove-close js-menu-close"></span>
          </div> <!-- /.select-menu-header -->

          <div class="select-menu-list js-navigation-container" role="menu">

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_included" name="do" type="radio" value="included" />
                <h4>Not watching</h4>
                <span class="description">You only receive notifications for conversations in which you participate or are @mentioned.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye-watch"></span>
                  Watch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item selected" role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input checked="checked" id="do_subscribed" name="do" type="radio" value="subscribed" />
                <h4>Watching</h4>
                <span class="description">You receive notifications for all conversations in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye-unwatch"></span>
                  Unwatch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_ignore" name="do" type="radio" value="ignore" />
                <h4>Ignoring</h4>
                <span class="description">You do not receive any notifications for conversations in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-mute"></span>
                  Stop ignoring
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

          </div> <!-- /.select-menu-list -->

        </div> <!-- /.select-menu-modal -->
      </div> <!-- /.select-menu-modal-holder -->
    </div> <!-- /.select-menu -->

</form>
    </li>

  <li>
  

  <div class="js-toggler-container js-social-container starring-container ">
    <a href="/ankit--sethi/GUI/unstar"
      class="minibutton with-count js-toggler-target star-button starred upwards"
      title="Unstar this repository" data-remote="true" data-method="post" rel="nofollow">
      <span class="octicon octicon-star-delete"></span><span class="text">Unstar</span>
    </a>

    <a href="/ankit--sethi/GUI/star"
      class="minibutton with-count js-toggler-target star-button unstarred upwards"
      title="Star this repository" data-remote="true" data-method="post" rel="nofollow">
      <span class="octicon octicon-star"></span><span class="text">Star</span>
    </a>

      <a class="social-count js-social-count" href="/ankit--sethi/GUI/stargazers">
        0
      </a>
  </div>

  </li>


        <li>
          <a href="/ankit--sethi/GUI/fork" class="minibutton with-count js-toggler-target fork-button lighter upwards" title="Fork this repo" rel="nofollow" data-method="post">
            <span class="octicon octicon-git-branch-create"></span><span class="text">Fork</span>
          </a>
          <a href="/ankit--sethi/GUI/network" class="social-count">51</a>
        </li>


</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo-forked"></span>
          <span class="author">
            <a href="/ankit--sethi" class="url fn" itemprop="url" rel="author"><span itemprop="title">ankit--sethi</span></a>
          </span>
          <span class="repohead-name-divider">/</span>
          <strong><a href="/ankit--sethi/GUI" class="js-current-repository js-repo-home-link">GUI</a></strong>

          <span class="page-context-loader">
            <img alt="Octocat-spinner-32" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

            <span class="fork-flag">
              <span class="text">forked from <a href="/open-ephys/GUI">open-ephys/GUI</a></span>
            </span>
        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">
      

      <div class="repository-with-sidebar repo-container new-discussion-timeline js-new-discussion-timeline  ">
        <div class="repository-sidebar">
            

<div class="sunken-menu vertical-right repo-nav js-repo-nav js-repository-container-pjax js-octicon-loaders">
  <div class="sunken-menu-contents">
    <ul class="sunken-menu-group">
      <li class="tooltipped leftwards" title="Code">
        <a href="/ankit--sethi/GUI" aria-label="Code" class="selected js-selected-navigation-item sunken-menu-item" data-gotokey="c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_tags repo_branches /ankit--sethi/GUI">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


      <li class="tooltipped leftwards" title="Pull Requests">
        <a href="/ankit--sethi/GUI/pulls" aria-label="Pull Requests" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-gotokey="p" data-selected-links="repo_pulls /ankit--sethi/GUI/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped leftwards" title="Wiki">
          <a href="/ankit--sethi/GUI/wiki" aria-label="Wiki" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_wiki /ankit--sethi/GUI/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="sunken-menu-separator"></div>
    <ul class="sunken-menu-group">

      <li class="tooltipped leftwards" title="Pulse">
        <a href="/ankit--sethi/GUI/pulse" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="pulse /ankit--sethi/GUI/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Graphs">
        <a href="/ankit--sethi/GUI/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_graphs repo_contributors /ankit--sethi/GUI/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped leftwards" title="Network">
        <a href="/ankit--sethi/GUI/network" aria-label="Network" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-selected-links="repo_network /ankit--sethi/GUI/network">
          <span class="octicon octicon-git-branch"></span> <span class="full-word">Network</span>
          <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
    </ul>


      <div class="sunken-menu-separator"></div>
      <ul class="sunken-menu-group">
        <li class="tooltipped leftwards" title="Settings">
          <a href="/ankit--sethi/GUI/settings"
            class="sunken-menu-item" data-pjax aria-label="Settings">
            <span class="octicon octicon-tools"></span> <span class="full-word">Settings</span>
            <img alt="Octocat-spinner-32" class="mini-loader" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32.gif" width="16" />
          </a>
        </li>
      </ul>
  </div>
</div>

              <div class="only-with-full-nav">
                

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=push">
  <h3><strong>HTTPS</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/ankit--sethi/GUI.git" readonly="readonly">

    <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/ankit--sethi/GUI.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="ssh"
  data-url="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=push">
  <h3><strong>SSH</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="git@github.com:ankit--sethi/GUI.git" readonly="readonly">

    <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="git@github.com:ankit--sethi/GUI.git" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=push">
  <h3><strong>Subversion</strong> checkout URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/ankit--sethi/GUI" readonly="readonly">

    <span class="js-zeroclipboard url-box-clippy minibutton zeroclipboard-button" data-clipboard-text="https://github.com/ankit--sethi/GUI" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


<p class="clone-options">You can clone with
      <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
      <a href="#" class="js-clone-selector" data-protocol="ssh">SSH</a>,
      or <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>.
  <span class="octicon help tooltipped upwards" title="Get help on which URL is right for you.">
    <a href="https://help.github.com/articles/which-remote-url-should-i-use">
    <span class="octicon octicon-question"></span>
    </a>
  </span>
</p>



                <a href="/ankit--sethi/GUI/archive/master.zip"
                   class="minibutton sidebar-button"
                   title="Download this repository as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
              </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<!-- blob contrib key: blob_contributors:v21:b672f8579650a91aabc71ac80a2a486f -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<a href="/ankit--sethi/GUI/find/master" data-pjax data-hotkey="t" class="js-show-file-finder" style="display:none">Show File Finder</a>

<div class="file-navigation">
  

<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target" data-hotkey="w"
    data-master-branch="master"
    data-ref="master"
    role="button" aria-label="Switch branches or tags" tabindex="0">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax>

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-remove-close js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Find or create a branch…" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Find or create a branch…">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/filterDirect/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="filterDirect"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="filterDirect">filterDirect</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/filterstuff/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="filterstuff"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="filterstuff">filterstuff</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/juce2/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="juce2"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="juce2">juce2</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/master/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="master"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="master">master</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/ripple/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="ripple"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="ripple">ripple</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/save_settings/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="save_settings"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="save_settings">save_settings</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/spikeplot/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="spikeplot"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="spikeplot">spikeplot</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/udpsource/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="udpsource"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="udpsource">udpsource</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/utilities/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="utilities"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="utilities">utilities</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/vs2012/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="vs2012"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="vs2012">vs2012</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/blob/windows/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="windows"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="windows">windows</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <form accept-charset="UTF-8" action="/ankit--sethi/GUI/branches" class="js-create-branch select-menu-item select-menu-new-item-form js-navigation-item js-new-item-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="wDXFYi/4cJYkfprIO/l9R4DNrLucl0avlTBH86JJUgQ=" /></div>
            <span class="octicon octicon-git-branch-create select-menu-item-icon"></span>
            <div class="select-menu-item-text">
              <h4>Create branch: <span class="js-new-item-name"></span></h4>
              <span class="description">from ‘master’</span>
            </div>
            <input type="hidden" name="name" id="name" class="js-new-item-value">
            <input type="hidden" name="branch" id="branch" value="master" />
            <input type="hidden" name="path" id="path" value="Source/Processors/Editors/ElectrodeButtons.cpp" />
          </form> <!-- /.select-menu-item -->

      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/tree/v0.1/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="v0.1"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="v0.1">v0.1</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/ankit--sethi/GUI/tree/v0.0/Source/Processors/Editors/ElectrodeButtons.cpp"
                 data-name="v0.0"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="v0.0">v0.0</a>
            </div> <!-- /.select-menu-item -->
        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/ankit--sethi/GUI" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">GUI</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/ankit--sethi/GUI/tree/master/Source" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">Source</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/ankit--sethi/GUI/tree/master/Source/Processors" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">Processors</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/ankit--sethi/GUI/tree/master/Source/Processors/Editors" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">Editors</span></a></span><span class="separator"> / </span><strong class="final-path">ElectrodeButtons.cpp</strong> <span class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="Source/Processors/Editors/ElectrodeButtons.cpp" data-copied-hint="copied!" title="copy to clipboard"><span class="octicon octicon-clippy"></span></span>
  </div>
</div>


  <div class="commit commit-loader file-history-tease js-deferred-content" data-url="/ankit--sethi/GUI/contributors/master/Source/Processors/Editors/ElectrodeButtons.cpp">
    Fetching contributors…

    <div class="participation">
      <p class="loader-loading"><img alt="Octocat-spinner-32-eaf2f5" height="16" src="https://github.global.ssl.fastly.net/images/spinners/octocat-spinner-32-EAF2F5.gif" width="16" /></p>
      <p class="loader-error">Cannot retrieve contributors at this time</p>
    </div>
  </div>

<div id="files" class="bubble">
  <div class="file">
    <div class="meta">
      <div class="info">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
          <span>166 lines (117 sloc)</span>
        <span>3.968 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
                <a class="minibutton js-update-url-with-hash"
                   href="/ankit--sethi/GUI/edit/master/Source/Processors/Editors/ElectrodeButtons.cpp"
                   data-method="post" rel="nofollow" data-hotkey="e">Edit</a>
          <a href="/ankit--sethi/GUI/raw/master/Source/Processors/Editors/ElectrodeButtons.cpp" class="button minibutton " id="raw-url">Raw</a>
            <a href="/ankit--sethi/GUI/blame/master/Source/Processors/Editors/ElectrodeButtons.cpp" class="button minibutton js-update-url-with-hash">Blame</a>
          <a href="/ankit--sethi/GUI/commits/master/Source/Processors/Editors/ElectrodeButtons.cpp" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->
          <a class="minibutton danger empty-icon tooltipped downwards"
             href="/ankit--sethi/GUI/delete/master/Source/Processors/Editors/ElectrodeButtons.cpp"
             title=""
             data-method="post" data-test-id="delete-blob-file" rel="nofollow">
          Delete
        </a>
      </div><!-- /.actions -->
    </div>
        <div class="blob-wrapper data type-c js-blob-data">
        <table class="file-code file-diff tab-size-8">
          <tr class="file-code-line">
            <td class="blob-line-nums">
              <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>

            </td>
            <td class="blob-line-code"><div class="code-body highlight"><pre><div class='line' id='LC1'><span class="cm">/*</span></div><div class='line' id='LC2'><span class="cm">    ------------------------------------------------------------------</span></div><div class='line' id='LC3'><br/></div><div class='line' id='LC4'><span class="cm">    This file is part of the Open Ephys GUI</span></div><div class='line' id='LC5'><span class="cm">    Copyright (C) 2013 Open Ephys</span></div><div class='line' id='LC6'><br/></div><div class='line' id='LC7'><span class="cm">    ------------------------------------------------------------------</span></div><div class='line' id='LC8'><br/></div><div class='line' id='LC9'><span class="cm">    This program is free software: you can redistribute it and/or modify</span></div><div class='line' id='LC10'><span class="cm">    it under the terms of the GNU General Public License as published by</span></div><div class='line' id='LC11'><span class="cm">    the Free Software Foundation, either version 3 of the License, or</span></div><div class='line' id='LC12'><span class="cm">    (at your option) any later version.</span></div><div class='line' id='LC13'><br/></div><div class='line' id='LC14'><span class="cm">    This program is distributed in the hope that it will be useful,</span></div><div class='line' id='LC15'><span class="cm">    but WITHOUT ANY WARRANTY; without even the implied warranty of</span></div><div class='line' id='LC16'><span class="cm">    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the</span></div><div class='line' id='LC17'><span class="cm">    GNU General Public License for more details.</span></div><div class='line' id='LC18'><br/></div><div class='line' id='LC19'><span class="cm">    You should have received a copy of the GNU General Public License</span></div><div class='line' id='LC20'><span class="cm">    along with this program.  If not, see &lt;http://www.gnu.org/licenses/&gt;.</span></div><div class='line' id='LC21'><br/></div><div class='line' id='LC22'><span class="cm">*/</span></div><div class='line' id='LC23'><br/></div><div class='line' id='LC24'><span class="cp">#include &quot;ElectrodeButtons.h&quot;</span></div><div class='line' id='LC25'><br/></div><div class='line' id='LC26'><span class="kt">void</span> <span class="n">ElectrodeButton</span><span class="o">::</span><span class="n">paintButton</span><span class="p">(</span><span class="n">Graphics</span><span class="o">&amp;</span> <span class="n">g</span><span class="p">,</span> <span class="kt">bool</span> <span class="n">isMouseOver</span><span class="p">,</span> <span class="kt">bool</span> <span class="n">isButtonDown</span><span class="p">)</span></div><div class='line' id='LC27'><span class="p">{</span></div><div class='line' id='LC28'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">getToggleState</span><span class="p">()</span> <span class="o">==</span> <span class="nb">true</span><span class="p">)</span></div><div class='line' id='LC29'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">orange</span><span class="p">);</span></div><div class='line' id='LC30'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span></div><div class='line' id='LC31'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">darkgrey</span><span class="p">);</span></div><div class='line' id='LC32'><br/></div><div class='line' id='LC33'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isMouseOver</span><span class="p">)</span></div><div class='line' id='LC34'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">white</span><span class="p">);</span></div><div class='line' id='LC35'><br/></div><div class='line' id='LC36'>	<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">isEnabled</span><span class="p">())</span></div><div class='line' id='LC37'>		<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">black</span><span class="p">);</span></div><div class='line' id='LC38'><br/></div><div class='line' id='LC39'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">fillRect</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">getWidth</span><span class="p">(),</span><span class="n">getHeight</span><span class="p">());</span></div><div class='line' id='LC40'><br/></div><div class='line' id='LC41'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// g.setFont(buttonFont);</span></div><div class='line' id='LC42'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">black</span><span class="p">);</span></div><div class='line' id='LC43'><br/></div><div class='line' id='LC44'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">drawRect</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">getWidth</span><span class="p">(),</span><span class="n">getHeight</span><span class="p">(),</span><span class="mf">1.0</span><span class="p">);</span></div><div class='line' id='LC45'><br/></div><div class='line' id='LC46'>	<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">isEnabled</span><span class="p">())</span></div><div class='line' id='LC47'>	<span class="p">{</span></div><div class='line' id='LC48'>		<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">grey</span><span class="p">);</span></div><div class='line' id='LC49'>	<span class="p">}</span></div><div class='line' id='LC50'><br/></div><div class='line' id='LC51'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">chan</span> <span class="o">&gt;=</span> <span class="mi">0</span><span class="p">)</span></div><div class='line' id='LC52'>		<span class="n">g</span><span class="p">.</span><span class="n">drawText</span><span class="p">(</span><span class="n">getButtonText</span><span class="p">(),</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">getWidth</span><span class="p">(),</span><span class="n">getHeight</span><span class="p">(),</span><span class="n">Justification</span><span class="o">::</span><span class="n">centred</span><span class="p">,</span><span class="nb">true</span><span class="p">);</span></div><div class='line' id='LC53'><span class="p">}</span></div><div class='line' id='LC54'><br/></div><div class='line' id='LC55'><span class="kt">void</span> <span class="n">ElectrodeButton</span><span class="o">::</span><span class="n">setChannelNum</span><span class="p">(</span><span class="kt">int</span> <span class="n">i</span><span class="p">)</span></div><div class='line' id='LC56'><span class="p">{</span></div><div class='line' id='LC57'>		<span class="n">setChannelNum</span><span class="p">(</span><span class="n">i</span><span class="p">,</span><span class="nb">true</span><span class="p">);</span></div><div class='line' id='LC58'><span class="p">}</span></div><div class='line' id='LC59'><br/></div><div class='line' id='LC60'><span class="kt">void</span> <span class="n">ElectrodeButton</span><span class="o">::</span><span class="n">setChannelNum</span><span class="p">(</span><span class="kt">int</span> <span class="n">i</span><span class="p">,</span> <span class="kt">bool</span> <span class="n">changeButtonText</span><span class="p">)</span></div><div class='line' id='LC61'><span class="p">{</span></div><div class='line' id='LC62'>	<span class="n">chan</span> <span class="o">=</span> <span class="n">i</span><span class="p">;</span></div><div class='line' id='LC63'><br/></div><div class='line' id='LC64'>	<span class="k">if</span> <span class="p">(</span><span class="n">changeButtonText</span><span class="p">)</span></div><div class='line' id='LC65'>	<span class="p">{</span></div><div class='line' id='LC66'>		<span class="n">setButtonText</span><span class="p">(</span><span class="n">String</span><span class="p">(</span><span class="n">chan</span><span class="p">));</span></div><div class='line' id='LC67'>	<span class="p">}</span></div><div class='line' id='LC68'><span class="p">}</span></div><div class='line' id='LC69'><br/></div><div class='line' id='LC70'><br/></div><div class='line' id='LC71'><span class="kt">void</span> <span class="n">ElectrodeEditorButton</span><span class="o">::</span><span class="n">paintButton</span><span class="p">(</span><span class="n">Graphics</span><span class="o">&amp;</span> <span class="n">g</span><span class="p">,</span> <span class="kt">bool</span> <span class="n">isMouseOver</span><span class="p">,</span> <span class="kt">bool</span> <span class="n">isButtonDown</span><span class="p">)</span></div><div class='line' id='LC72'><span class="p">{</span></div><div class='line' id='LC73'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">getToggleState</span><span class="p">()</span> <span class="o">==</span> <span class="nb">true</span><span class="p">)</span></div><div class='line' id='LC74'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">darkgrey</span><span class="p">);</span></div><div class='line' id='LC75'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span></div><div class='line' id='LC76'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">lightgrey</span><span class="p">);</span></div><div class='line' id='LC77'><br/></div><div class='line' id='LC78'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setFont</span><span class="p">(</span><span class="n">font</span><span class="p">);</span></div><div class='line' id='LC79'><br/></div><div class='line' id='LC80'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">drawText</span><span class="p">(</span><span class="n">name</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="mi">0</span><span class="p">,</span><span class="n">getWidth</span><span class="p">(),</span><span class="n">getHeight</span><span class="p">(),</span><span class="n">Justification</span><span class="o">::</span><span class="n">left</span><span class="p">,</span><span class="nb">true</span><span class="p">);</span></div><div class='line' id='LC81'><span class="p">}</span></div><div class='line' id='LC82'><br/></div><div class='line' id='LC83'><span class="n">ThresholdSlider</span><span class="o">::</span><span class="n">ThresholdSlider</span><span class="p">(</span><span class="n">Font</span> <span class="n">f</span><span class="p">)</span> <span class="o">:</span> <span class="n">Slider</span><span class="p">(</span><span class="s">&quot;name&quot;</span><span class="p">),</span> <span class="n">font</span><span class="p">(</span><span class="n">f</span><span class="p">)</span></div><div class='line' id='LC84'><span class="p">{</span></div><div class='line' id='LC85'><br/></div><div class='line' id='LC86'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">setSliderStyle</span><span class="p">(</span><span class="n">Slider</span><span class="o">::</span><span class="n">Rotary</span><span class="p">);</span></div><div class='line' id='LC87'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">setRange</span><span class="p">(</span><span class="mf">25.0f</span><span class="p">,</span><span class="mf">400.0f</span><span class="p">,</span><span class="mf">25.0f</span><span class="p">);</span></div><div class='line' id='LC88'>&nbsp;&nbsp;&nbsp;<span class="c1">// setValue(75.0f);</span></div><div class='line' id='LC89'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">setTextBoxStyle</span><span class="p">(</span><span class="n">Slider</span><span class="o">::</span><span class="n">NoTextBox</span><span class="p">,</span> <span class="nb">false</span><span class="p">,</span> <span class="mi">40</span><span class="p">,</span> <span class="mi">20</span><span class="p">);</span></div><div class='line' id='LC90'><br/></div><div class='line' id='LC91'><span class="p">}</span></div><div class='line' id='LC92'><br/></div><div class='line' id='LC93'><span class="kt">void</span> <span class="n">ThresholdSlider</span><span class="o">::</span><span class="n">paint</span><span class="p">(</span><span class="n">Graphics</span><span class="o">&amp;</span> <span class="n">g</span><span class="p">)</span></div><div class='line' id='LC94'><span class="p">{</span></div><div class='line' id='LC95'><br/></div><div class='line' id='LC96'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ColourGradient</span> <span class="n">grad</span> <span class="o">=</span> <span class="n">ColourGradient</span><span class="p">(</span><span class="n">Colour</span><span class="p">(</span><span class="mi">40</span><span class="p">,</span> <span class="mi">40</span><span class="p">,</span> <span class="mi">40</span><span class="p">),</span> <span class="mf">0.0f</span><span class="p">,</span> <span class="mf">0.0f</span><span class="p">,</span></div><div class='line' id='LC97'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Colour</span><span class="p">(</span><span class="mi">80</span><span class="p">,</span> <span class="mi">80</span><span class="p">,</span> <span class="mi">80</span><span class="p">),</span> <span class="mf">0.0</span><span class="p">,</span> <span class="mf">40.0f</span><span class="p">,</span> <span class="nb">false</span><span class="p">);</span></div><div class='line' id='LC98'><br/></div><div class='line' id='LC99'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Path</span> <span class="n">p</span><span class="p">;</span></div><div class='line' id='LC100'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">p</span><span class="p">.</span><span class="n">addPieSegment</span><span class="p">(</span><span class="mi">3</span><span class="p">,</span> <span class="mi">3</span><span class="p">,</span> <span class="n">getWidth</span><span class="p">()</span><span class="o">-</span><span class="mi">6</span><span class="p">,</span> <span class="n">getHeight</span><span class="p">()</span><span class="o">-</span><span class="mi">6</span><span class="p">,</span> <span class="mi">5</span><span class="o">*</span><span class="n">double_Pi</span><span class="o">/</span><span class="mi">4</span><span class="o">-</span><span class="mf">0.2</span><span class="p">,</span> <span class="mi">5</span><span class="o">*</span><span class="n">double_Pi</span><span class="o">/</span><span class="mi">4</span><span class="o">+</span><span class="mi">3</span><span class="o">*</span><span class="n">double_Pi</span><span class="o">/</span><span class="mi">2</span><span class="o">+</span><span class="mf">0.2</span><span class="p">,</span> <span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC101'><br/></div><div class='line' id='LC102'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setGradientFill</span><span class="p">(</span><span class="n">grad</span><span class="p">);</span></div><div class='line' id='LC103'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">fillPath</span><span class="p">(</span><span class="n">p</span><span class="p">);</span></div><div class='line' id='LC104'><br/></div><div class='line' id='LC105'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">String</span> <span class="n">valueString</span><span class="p">;</span></div><div class='line' id='LC106'><br/></div><div class='line' id='LC107'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">isActive</span><span class="p">)</span></div><div class='line' id='LC108'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC109'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">p</span> <span class="o">=</span> <span class="n">makeRotaryPath</span><span class="p">(</span><span class="n">getMinimum</span><span class="p">(),</span> <span class="n">getMaximum</span><span class="p">(),</span> <span class="n">getValue</span><span class="p">());</span></div><div class='line' id='LC110'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colour</span><span class="p">(</span><span class="mi">240</span><span class="p">,</span><span class="mi">179</span><span class="p">,</span><span class="mi">12</span><span class="p">));</span></div><div class='line' id='LC111'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">fillPath</span><span class="p">(</span><span class="n">p</span><span class="p">);</span></div><div class='line' id='LC112'><br/></div><div class='line' id='LC113'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">valueString</span> <span class="o">=</span> <span class="n">String</span><span class="p">((</span><span class="kt">int</span><span class="p">)</span> <span class="n">getValue</span><span class="p">());</span></div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC115'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span></div><div class='line' id='LC116'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC117'><br/></div><div class='line' id='LC118'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">valueString</span> <span class="o">=</span> <span class="s">&quot;&quot;</span><span class="p">;</span></div><div class='line' id='LC119'><br/></div><div class='line' id='LC120'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">valueArray</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC121'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC122'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">p</span> <span class="o">=</span> <span class="n">makeRotaryPath</span><span class="p">(</span><span class="n">getMinimum</span><span class="p">(),</span> <span class="n">getMaximum</span><span class="p">(),</span> <span class="n">valueArray</span><span class="p">[</span><span class="n">i</span><span class="p">]);</span></div><div class='line' id='LC123'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">lightgrey</span><span class="p">.</span><span class="n">withAlpha</span><span class="p">(</span><span class="mf">0.4f</span><span class="p">));</span></div><div class='line' id='LC124'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">fillPath</span><span class="p">(</span><span class="n">p</span><span class="p">);</span></div><div class='line' id='LC125'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">valueString</span> <span class="o">=</span> <span class="n">String</span><span class="p">((</span><span class="kt">int</span><span class="p">)</span> <span class="n">valueArray</span><span class="p">.</span><span class="n">getLast</span><span class="p">());</span></div><div class='line' id='LC126'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC127'><br/></div><div class='line' id='LC128'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC129'><br/></div><div class='line' id='LC130'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">font</span><span class="p">.</span><span class="n">setHeight</span><span class="p">(</span><span class="mf">9.0</span><span class="p">);</span></div><div class='line' id='LC131'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setFont</span><span class="p">(</span><span class="n">font</span><span class="p">);</span></div><div class='line' id='LC132'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">int</span> <span class="n">stringWidth</span> <span class="o">=</span> <span class="n">font</span><span class="p">.</span><span class="n">getStringWidth</span><span class="p">(</span><span class="n">valueString</span><span class="p">);</span></div><div class='line' id='LC133'><br/></div><div class='line' id='LC134'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setFont</span><span class="p">(</span><span class="n">font</span><span class="p">);</span></div><div class='line' id='LC135'><br/></div><div class='line' id='LC136'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">setColour</span><span class="p">(</span><span class="n">Colours</span><span class="o">::</span><span class="n">darkgrey</span><span class="p">);</span></div><div class='line' id='LC137'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">g</span><span class="p">.</span><span class="n">drawSingleLineText</span><span class="p">(</span><span class="n">valueString</span><span class="p">,</span> <span class="n">getWidth</span><span class="p">()</span><span class="o">/</span><span class="mi">2</span> <span class="o">-</span> <span class="n">stringWidth</span><span class="o">/</span><span class="mi">2</span><span class="p">,</span> <span class="n">getHeight</span><span class="p">()</span><span class="o">/</span><span class="mi">2</span><span class="o">+</span><span class="mi">3</span><span class="p">);</span></div><div class='line' id='LC138'><br/></div><div class='line' id='LC139'><span class="p">}</span></div><div class='line' id='LC140'><br/></div><div class='line' id='LC141'><span class="n">Path</span> <span class="n">ThresholdSlider</span><span class="o">::</span><span class="n">makeRotaryPath</span><span class="p">(</span><span class="kt">double</span> <span class="n">min</span><span class="p">,</span> <span class="kt">double</span> <span class="n">max</span><span class="p">,</span> <span class="kt">double</span> <span class="n">val</span><span class="p">)</span></div><div class='line' id='LC142'><span class="p">{</span></div><div class='line' id='LC143'><br/></div><div class='line' id='LC144'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Path</span> <span class="n">p</span><span class="p">;</span></div><div class='line' id='LC145'><br/></div><div class='line' id='LC146'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">start</span> <span class="o">=</span> <span class="mi">5</span><span class="o">*</span><span class="n">double_Pi</span><span class="o">/</span><span class="mi">4</span> <span class="o">-</span> <span class="mf">0.11</span><span class="p">;</span></div><div class='line' id='LC147'><br/></div><div class='line' id='LC148'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">range</span> <span class="o">=</span> <span class="p">(</span><span class="n">val</span><span class="o">-</span><span class="n">min</span><span class="p">)</span><span class="o">/</span><span class="p">(</span><span class="n">max</span> <span class="o">-</span> <span class="n">min</span><span class="p">)</span><span class="o">*</span><span class="mf">1.5</span><span class="o">*</span><span class="n">double_Pi</span> <span class="o">+</span> <span class="n">start</span> <span class="o">+</span> <span class="mf">0.22</span><span class="p">;</span></div><div class='line' id='LC149'><br/></div><div class='line' id='LC150'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">p</span><span class="p">.</span><span class="n">addPieSegment</span><span class="p">(</span><span class="mi">6</span><span class="p">,</span><span class="mi">6</span><span class="p">,</span> <span class="n">getWidth</span><span class="p">()</span><span class="o">-</span><span class="mi">12</span><span class="p">,</span> <span class="n">getHeight</span><span class="p">()</span><span class="o">-</span><span class="mi">12</span><span class="p">,</span> <span class="n">start</span><span class="p">,</span> <span class="n">range</span><span class="p">,</span> <span class="mf">0.65</span><span class="p">);</span></div><div class='line' id='LC151'><br/></div><div class='line' id='LC152'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">return</span> <span class="n">p</span><span class="p">;</span></div><div class='line' id='LC153'><br/></div><div class='line' id='LC154'><span class="p">}</span></div><div class='line' id='LC155'><br/></div><div class='line' id='LC156'><span class="kt">void</span> <span class="n">ThresholdSlider</span><span class="o">::</span><span class="n">setActive</span><span class="p">(</span><span class="kt">bool</span> <span class="n">t</span><span class="p">)</span></div><div class='line' id='LC157'><span class="p">{</span></div><div class='line' id='LC158'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">isActive</span> <span class="o">=</span> <span class="n">t</span><span class="p">;</span></div><div class='line' id='LC159'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">repaint</span><span class="p">();</span></div><div class='line' id='LC160'><span class="p">}</span></div><div class='line' id='LC161'><br/></div><div class='line' id='LC162'><span class="kt">void</span> <span class="n">ThresholdSlider</span><span class="o">::</span><span class="n">setValues</span><span class="p">(</span><span class="n">Array</span><span class="o">&lt;</span><span class="kt">double</span><span class="o">&gt;</span> <span class="n">v</span><span class="p">)</span></div><div class='line' id='LC163'><span class="p">{</span></div><div class='line' id='LC164'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">valueArray</span> <span class="o">=</span> <span class="n">v</span><span class="p">;</span></div><div class='line' id='LC165'><span class="p">}</span></div></pre></div></td>
          </tr>
        </table>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github" title="GitHub"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2014 <span title="0.07714s from github-fe122-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
          <div class="suggester-container">
              <div class="suggester fullscreen-suggester js-navigation-container" id="fullscreen_suggester"
                 data-url="/ankit--sethi/GUI/suggestions/commit">
              </div>
          </div>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped leftwards" title="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped leftwards"
      title="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-remove-close close js-ajax-error-dismiss"></a>
      Something went wrong with that request. Please try again.
    </div>

  </body>
</html>

