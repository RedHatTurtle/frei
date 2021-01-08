<!DOCTYPE html>
<html class="" lang="en">
<head prefix="og: http://ogp.me/ns#">
<meta charset="utf-8">
<link as="style" href="https://gitlab.onelab.info/assets/application-051048a171ccf14f73419f46d3bd8204aa3ed585a72924faea0192f53d42cfce.css" rel="preload">
<link as="style" href="https://gitlab.onelab.info/assets/highlight/themes/white-aa4568025f9b4ea36b357bdccb95c9138a515f1e611b59f20a1777a68b6995db.css" rel="preload">
<link as="font" href="https://gitlab.onelab.info/assets/fontawesome-webfont-2adefcbc041e7d18fcf2d417879dc5a09997aa64d675b7a3c4b6ce33da13f3fe.woff2?v=4.7.0" rel="preload" type="font/woff2">

<meta content="IE=edge" http-equiv="X-UA-Compatible">

<meta content="object" property="og:type">
<meta content="GitLab" property="og:site_name">
<meta content="Common/GmshDefines.h · master · gmsh / gmsh" property="og:title">
<meta content="The main Gmsh project (http://gmsh.info)" property="og:description">
<meta content="https://gitlab.onelab.info/uploads/-/system/project/avatar/3/gmsh.png" property="og:image">
<meta content="64" property="og:image:width">
<meta content="64" property="og:image:height">
<meta content="https://gitlab.onelab.info/gmsh/gmsh/blob/master/Common/GmshDefines.h" property="og:url">
<meta content="summary" property="twitter:card">
<meta content="Common/GmshDefines.h · master · gmsh / gmsh" property="twitter:title">
<meta content="The main Gmsh project (http://gmsh.info)" property="twitter:description">
<meta content="https://gitlab.onelab.info/uploads/-/system/project/avatar/3/gmsh.png" property="twitter:image">

<title>Common/GmshDefines.h · master · gmsh / gmsh · GitLab</title>
<meta content="The main Gmsh project (http://gmsh.info)" name="description">
<link rel="shortcut icon" type="image/png" href="/assets/favicon-7901bd695fb93edb07975966062049829afb56cf11511236e61bcf425070e36e.png" id="favicon" data-original-href="/assets/favicon-7901bd695fb93edb07975966062049829afb56cf11511236e61bcf425070e36e.png" />

<link rel="stylesheet" media="all" href="/assets/application-051048a171ccf14f73419f46d3bd8204aa3ed585a72924faea0192f53d42cfce.css" />

<link rel="stylesheet" media="all" href="/assets/application_utilities-f636133344c99e44738666d9eecff2eb037c36ed05c5cbc2bfacc8dfad2bede2.css" />
<link rel="stylesheet" media="all" href="/assets/themes/theme_indigo-de5d6083627877d91d7026811b987cc8ef64f626958073109b5dfb8cb510a7b4.css" />

<link rel="stylesheet" media="all" href="/assets/highlight/themes/white-aa4568025f9b4ea36b357bdccb95c9138a515f1e611b59f20a1777a68b6995db.css" />


<script>
//<![CDATA[
window.gon={};gon.features={"suggestPipeline":true,"gitlabCiYmlPreview":false};
//]]>
</script>



<script src="/assets/webpack/runtime.95010ba7.bundle.js" defer="defer"></script>
<script src="/assets/webpack/main.88aefd49.chunk.js" defer="defer"></script>
<script src="/assets/webpack/commons-globalSearch-pages.admin.abuse_reports-pages.admin.groups.show-pages.admin.projects-pages.ad-b36b8820.a8bba176.chunk.js" defer="defer"></script>
<script src="/assets/webpack/commons-pages.groups.boards-pages.groups.details-pages.groups.show-pages.projects-pages.projects.act-d59c63d5.c5897a66.chunk.js" defer="defer"></script>
<script src="/assets/webpack/commons-pages.admin.application_settings-pages.admin.application_settings.general-pages.admin.applic-e5cd0c99.3135386d.chunk.js" defer="defer"></script>
<script src="/assets/webpack/commons-pages.groups.milestones.edit-pages.groups.milestones.new-pages.projects.blame.show-pages.pro-77e8c306.37c1e60d.chunk.js" defer="defer"></script>
<script src="/assets/webpack/commons-pages.projects.blob.show-pages.projects.shared.web_ide_link-pages.projects.show-pages.projec-6b390939.d5832f7d.chunk.js" defer="defer"></script>
<script src="/assets/webpack/commons-pages.projects.blame.show-pages.projects.blob.show.0a845ac7.chunk.js" defer="defer"></script>
<script src="/assets/webpack/pages.projects.blob.show.fd89f639.chunk.js" defer="defer"></script>


<meta name="csrf-param" content="authenticity_token" />
<meta name="csrf-token" content="cYVIZVm0qDtzraQ6eMFzam0dv3IVeaLYLZ6VLDm9jR0RTwnyzv+j+jMgMNcGQC3T39uuZ1/zDMKGcEamlys88Q==" />
<meta name="csp-nonce" />
<meta name="action-cable-url" content="/-/cable" />
<meta content="width=device-width, initial-scale=1, maximum-scale=1" name="viewport">
<meta content="#474D57" name="theme-color">
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-iphone-5a9cee0e8a51212e70b90c87c12f382c428870c0ff67d1eb034d884b78d2dae7.png" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-ipad-a6eec6aeb9da138e507593b464fdac213047e49d3093fc30e90d9a995df83ba3.png" sizes="76x76" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-iphone-retina-72e2aadf86513a56e050e7f0f2355deaa19cc17ed97bbe5147847f2748e5a3e3.png" sizes="120x120" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-ipad-retina-8ebe416f5313483d9c1bc772b5bbe03ecad52a54eba443e5215a22caed2a16a2.png" sizes="152x152" />
<link color="rgb(226, 67, 41)" href="/assets/logo-d36b5212042cebc89b96df4bf6ac24e43db316143e89926c0db839ff694d2de4.svg" rel="mask-icon">
<meta content="/assets/msapplication-tile-1196ec67452f618d39cdd85e2e3a542f76574c071051ae7effbfde01710eb17d.png" name="msapplication-TileImage">
<meta content="#30353E" name="msapplication-TileColor">




</head>

<body class="ui-indigo tab-width-8  gl-browser-generic gl-platform-other" data-find-file="/gmsh/gmsh/-/find_file/master" data-group="gmsh" data-namespace-id="2" data-page="projects:blob:show" data-page-type-id="master/Common/GmshDefines.h" data-project="gmsh" data-project-id="3">

<script>
//<![CDATA[
gl = window.gl || {};
gl.client = {"isGeneric":true,"isOther":true};


//]]>
</script>


<header class="navbar navbar-gitlab navbar-expand-sm js-navbar" data-qa-selector="navbar">
<a class="sr-only gl-accessibility" href="#content-body" tabindex="1">Skip to content</a>
<div class="container-fluid">
<div class="header-content">
<div class="title-container">
<h1 class="title">
<span class="gl-sr-only">GitLab</span>
<a title="Dashboard" id="logo" href="/"><svg width="24" height="24" class="tanuki-logo" viewBox="0 0 36 36">
  <path class="tanuki-shape tanuki-left-ear" fill="#e24329" d="M2 14l9.38 9v-9l-4-12.28c-.205-.632-1.176-.632-1.38 0z"/>
  <path class="tanuki-shape tanuki-right-ear" fill="#e24329" d="M34 14l-9.38 9v-9l4-12.28c.205-.632 1.176-.632 1.38 0z"/>
  <path class="tanuki-shape tanuki-nose" fill="#e24329" d="M18,34.38 3,14 33,14 Z"/>
  <path class="tanuki-shape tanuki-left-eye" fill="#fc6d26" d="M18,34.38 11.38,14 2,14 6,25Z"/>
  <path class="tanuki-shape tanuki-right-eye" fill="#fc6d26" d="M18,34.38 24.62,14 34,14 30,25Z"/>
  <path class="tanuki-shape tanuki-left-cheek" fill="#fca326" d="M2 14L.1 20.16c-.18.565 0 1.2.5 1.56l17.42 12.66z"/>
  <path class="tanuki-shape tanuki-right-cheek" fill="#fca326" d="M34 14l1.9 6.16c.18.565 0 1.2-.5 1.56L18 34.38z"/>
</svg>

<span class="logo-text d-none d-lg-block gl-ml-3">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 617 169"><path d="M315.26 2.97h-21.8l.1 162.5h88.3v-20.1h-66.5l-.1-142.4M465.89 136.95c-5.5 5.7-14.6 11.4-27 11.4-16.6 0-23.3-8.2-23.3-18.9 0-16.1 11.2-23.8 35-23.8 4.5 0 11.7.5 15.4 1.2v30.1h-.1m-22.6-98.5c-17.6 0-33.8 6.2-46.4 16.7l7.7 13.4c8.9-5.2 19.8-10.4 35.5-10.4 17.9 0 25.8 9.2 25.8 24.6v7.9c-3.5-.7-10.7-1.2-15.1-1.2-38.2 0-57.6 13.4-57.6 41.4 0 25.1 15.4 37.7 38.7 37.7 15.7 0 30.8-7.2 36-18.9l4 15.9h15.4v-83.2c-.1-26.3-11.5-43.9-44-43.9M557.63 149.1c-8.2 0-15.4-1-20.8-3.5V70.5c7.4-6.2 16.6-10.7 28.3-10.7 21.1 0 29.2 14.9 29.2 39 0 34.2-13.1 50.3-36.7 50.3m9.2-110.6c-19.5 0-30 13.3-30 13.3v-21l-.1-27.8h-21.3l.1 158.5c10.7 4.5 25.3 6.9 41.2 6.9 40.7 0 60.3-26 60.3-70.9-.1-35.5-18.2-59-50.2-59M77.9 20.6c19.3 0 31.8 6.4 39.9 12.9l9.4-16.3C114.5 6 97.3 0 78.9 0 32.5 0 0 28.3 0 85.4c0 59.8 35.1 83.1 75.2 83.1 20.1 0 37.2-4.7 48.4-9.4l-.5-63.9V75.1H63.6v20.1h38l.5 48.5c-5 2.5-13.6 4.5-25.3 4.5-32.2 0-53.8-20.3-53.8-63-.1-43.5 22.2-64.6 54.9-64.6M231.43 2.95h-21.3l.1 27.3v94.3c0 26.3 11.4 43.9 43.9 43.9 4.5 0 8.9-.4 13.1-1.2v-19.1c-3.1.5-6.4.7-9.9.7-17.9 0-25.8-9.2-25.8-24.6v-65h35.7v-17.8h-35.7l-.1-38.5M155.96 165.47h21.3v-124h-21.3v124M155.96 24.37h21.3V3.07h-21.3v21.3"/></svg>

</span>
</a></h1>
<ul class="list-unstyled navbar-sub-nav">
<li class="home"><a title="Projects" class="dashboard-shortcuts-projects" href="/explore">Projects
</a></li><li class=""><a title="Groups" class="dashboard-shortcuts-groups" href="/explore/groups">Groups
</a></li><li class=""><a title="Snippets" class="dashboard-shortcuts-snippets" href="/explore/snippets">Snippets
</a></li><li>
<a title="About GitLab CE" href="/help">Help</a>
</li>
</ul>

</div>
<div class="navbar-collapse collapse">
<ul class="nav navbar-nav">
<li class="nav-item d-none d-lg-block m-auto">
<div class="search search-form" data-track-event="activate_form_input" data-track-label="navbar_search" data-track-value="">
<form class="form-inline" action="/search" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" /><div class="search-input-container">
<div class="search-input-wrap">
<div class="dropdown" data-url="/search/autocomplete">
<input type="search" name="search" id="search" placeholder="Search or jump to…" class="search-input dropdown-menu-toggle no-outline js-search-dashboard-options" spellcheck="false" autocomplete="off" data-issues-path="/dashboard/issues" data-mr-path="/dashboard/merge_requests" data-qa-selector="search_term_field" aria-label="Search or jump to…" />
<button class="hidden js-dropdown-search-toggle" data-toggle="dropdown" type="button"></button>
<div class="dropdown-menu dropdown-select" data-testid="dashboard-search-options">
<div class="dropdown-content"><ul>
<li class="dropdown-menu-empty-item">
<a>
Loading...
</a>
</li>
</ul>
</div><div class="dropdown-loading"><div class="gl-spinner-container"><span class="gl-spinner gl-spinner-orange gl-spinner-md gl-mt-7" aria-label="Loading"></span></div></div>
</div>
<svg class="s16 search-icon" data-testid="search-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#search"></use></svg>
<svg class="s16 clear-icon js-clear-input" data-testid="close-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#close"></use></svg>
</div>
</div>
</div>
<input type="hidden" name="group_id" id="group_id" value="2" class="js-search-group-options" data-group-path="gmsh" data-name="gmsh" data-issues-path="/groups/gmsh/-/issues" data-mr-path="/groups/gmsh/-/merge_requests" />
<input type="hidden" name="project_id" id="search_project_id" value="3" class="js-search-project-options" data-project-path="gmsh" data-name="gmsh" data-issues-path="/gmsh/gmsh/-/issues" data-mr-path="/gmsh/gmsh/-/merge_requests" data-issues-disabled="false" />
<input type="hidden" name="scope" id="scope" />
<input type="hidden" name="search_code" id="search_code" value="true" />
<input type="hidden" name="snippets" id="snippets" value="false" />
<input type="hidden" name="repository_ref" id="repository_ref" value="master" />
<input type="hidden" name="nav_source" id="nav_source" value="navbar" />
<div class="search-autocomplete-opts hide" data-autocomplete-path="/search/autocomplete" data-autocomplete-project-id="3" data-autocomplete-project-ref="master"></div>
</form></div>

</li>
<li class="nav-item d-inline-block d-lg-none">
<a title="Search" aria-label="Search" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/search?project_id=3"><svg class="s16" data-testid="search-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#search"></use></svg>
</a></li>
<li class="nav-item header-help dropdown d-none d-md-block">
<a class="header-help-dropdown-toggle" data-toggle="dropdown" href="/help"><span class="gl-sr-only">
Help
</span>
<svg class="s16" data-testid="question-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#question"></use></svg>
<svg class="s16 caret-down" data-testid="chevron-down-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#chevron-down"></use></svg>
</a><div class="dropdown-menu dropdown-menu-right">
<ul>

<li>
<a href="/help">Help</a>
</li>
<li>
<a href="https://about.gitlab.com/getting-help/">Support</a>
</li>
<li>
<a target="_blank" class="text-nowrap" rel="noopener noreferrer" data-track-event="click_forum" data-track-property="question_menu" href="https://forum.gitlab.com/">Community forum</a>

</li>
<li>
<button class="js-shortcuts-modal-trigger" type="button">
Keyboard shortcuts
<span aria-hidden class="text-secondary float-right">?</span>
</button>
</li>
<li class="divider"></li>
<li>
<a href="https://about.gitlab.com/submit-feedback">Submit feedback</a>
</li>
<li>
<a target="_blank" class="text-nowrap" href="https://about.gitlab.com/contributing">Contribute to GitLab
</a>
</li>

</ul>

</div>
</li>
<li class="nav-item">
<div>
<a class="btn btn-sign-in" href="/users/sign_in?redirect_to_referer=yes">Sign in / Register</a>
</div>
</li>
</ul>
</div>
<button class="navbar-toggler d-block d-sm-none" type="button">
<span class="sr-only">Toggle navigation</span>
<svg class="s12 more-icon js-navbar-toggle-right" data-testid="ellipsis_h-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#ellipsis_h"></use></svg>
<svg class="s12 close-icon js-navbar-toggle-left" data-testid="close-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#close"></use></svg>
</button>
</div>
</div>
</header>

<div class="layout-page page-with-contextual-sidebar">
<div class="nav-sidebar">
<div class="nav-sidebar-inner-scroll">
<div class="context-header">
<a title="gmsh" href="/gmsh/gmsh"><div class="avatar-container rect-avatar s40 project-avatar">
<img alt="gmsh" class="avatar s40 avatar-tile lazy" width="40" height="40" data-src="/uploads/-/system/project/avatar/3/gmsh.png" src="data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEKAAEALAAAAAABAAEAAAICTAEAOw==" />
</div>
<div class="sidebar-context-title">
gmsh
</div>
</a></div>
<ul class="sidebar-top-level-items qa-project-sidebar">
<li class="home"><a class="shortcuts-project rspec-project-link" data-qa-selector="project_link" href="/gmsh/gmsh"><div class="nav-icon-container">
<svg class="s16" data-testid="home-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#home"></use></svg>
</div>
<span class="nav-item-name">
Project overview
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item"><a href="/gmsh/gmsh"><strong class="fly-out-top-item-name">
Project overview
</strong>
</a></li><li class="divider fly-out-top-item"></li>
<li class=""><a title="Project details" class="shortcuts-project" href="/gmsh/gmsh"><span>Details</span>
</a></li><li class=""><a title="Activity" class="shortcuts-project-activity" data-qa-selector="activity_link" href="/gmsh/gmsh/activity"><span>Activity</span>
</a></li><li class=""><a title="Releases" class="shortcuts-project-releases" href="/gmsh/gmsh/-/releases"><span>Releases</span>
</a></li></ul>
</li><li class="active"><a class="shortcuts-tree" data-qa-selector="repository_link" href="/gmsh/gmsh/-/tree/master"><div class="nav-icon-container">
<svg class="s16" data-testid="doc-text-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#doc-text"></use></svg>
</div>
<span class="nav-item-name" id="js-onboarding-repo-link">
Repository
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item active"><a href="/gmsh/gmsh/-/tree/master"><strong class="fly-out-top-item-name">
Repository
</strong>
</a></li><li class="divider fly-out-top-item"></li>
<li class="active"><a href="/gmsh/gmsh/-/tree/master">Files
</a></li><li class=""><a id="js-onboarding-commits-link" href="/gmsh/gmsh/-/commits/master">Commits
</a></li><li class=""><a data-qa-selector="branches_link" id="js-onboarding-branches-link" href="/gmsh/gmsh/-/branches">Branches
</a></li><li class=""><a data-qa-selector="tags_link" href="/gmsh/gmsh/-/tags">Tags
</a></li><li class=""><a href="/gmsh/gmsh/-/graphs/master">Contributors
</a></li><li class=""><a href="/gmsh/gmsh/-/network/master">Graph
</a></li><li class=""><a href="/gmsh/gmsh/-/compare?from=master&amp;to=master">Compare
</a></li>
</ul>
</li><li class=""><a class="shortcuts-issues qa-issues-item" href="/gmsh/gmsh/-/issues"><div class="nav-icon-container">
<svg class="s16" data-testid="issues-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#issues"></use></svg>
</div>
<span class="nav-item-name" id="js-onboarding-issues-link">
Issues
</span>
<span class="badge badge-pill count issue_counter">
106
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item"><a href="/gmsh/gmsh/-/issues"><strong class="fly-out-top-item-name">
Issues
</strong>
<span class="badge badge-pill count issue_counter fly-out-badge">
106
</span>
</a></li><li class="divider fly-out-top-item"></li>
<li class=""><a title="Issues" href="/gmsh/gmsh/-/issues"><span>
List
</span>
</a></li><li class=""><a title="Boards" data-qa-selector="issue_boards_link" href="/gmsh/gmsh/-/boards"><span>
Boards
</span>
</a></li><li class=""><a title="Labels" class="qa-labels-link" href="/gmsh/gmsh/-/labels"><span>
Labels
</span>
</a></li><li class=""><a title="Service Desk" href="/gmsh/gmsh/-/issues/service_desk">Service Desk
</a></li>
<li class=""><a title="Milestones" class="qa-milestones-link" href="/gmsh/gmsh/-/milestones"><span>
Milestones
</span>
</a></li>
</ul>
</li><li class=""><a class="shortcuts-merge_requests" data-qa-selector="merge_requests_link" href="/gmsh/gmsh/-/merge_requests"><div class="nav-icon-container">
<svg class="s16" data-testid="git-merge-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#git-merge"></use></svg>
</div>
<span class="nav-item-name" id="js-onboarding-mr-link">
Merge Requests
</span>
<span class="badge badge-pill count merge_counter js-merge-counter">
7
</span>
</a><ul class="sidebar-sub-level-items is-fly-out-only">
<li class="fly-out-top-item"><a href="/gmsh/gmsh/-/merge_requests"><strong class="fly-out-top-item-name">
Merge Requests
</strong>
<span class="badge badge-pill count merge_counter js-merge-counter fly-out-badge">
7
</span>
</a></li></ul>
</li>
<li class=""><a class="shortcuts-pipelines qa-link-pipelines rspec-link-pipelines" data-qa-selector="ci_cd_link" href="/gmsh/gmsh/-/pipelines"><div class="nav-icon-container">
<svg class="s16" data-testid="rocket-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#rocket"></use></svg>
</div>
<span class="nav-item-name" id="js-onboarding-pipelines-link">
CI / CD
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item"><a href="/gmsh/gmsh/-/pipelines"><strong class="fly-out-top-item-name">
CI / CD
</strong>
</a></li><li class="divider fly-out-top-item"></li>
<li class=""><a title="Pipelines" class="shortcuts-pipelines" href="/gmsh/gmsh/-/pipelines"><span>
Pipelines
</span>
</a></li><li class=""><a title="Jobs" class="shortcuts-builds" href="/gmsh/gmsh/-/jobs"><span>
Jobs
</span>
</a></li><li class=""><a title="Schedules" class="shortcuts-builds" href="/gmsh/gmsh/-/pipeline_schedules"><span>
Schedules
</span>
</a></li>
</ul>
</li>
<li class=""><a class="shortcuts-operations" data-qa-selector="operations_link" href="/gmsh/gmsh/-/environments/metrics"><div class="nav-icon-container">
<svg class="s16" data-testid="cloud-gear-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#cloud-gear"></use></svg>
</div>
<span class="nav-item-name">
Operations
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item"><a href="/gmsh/gmsh/-/environments/metrics"><strong class="fly-out-top-item-name">
Operations
</strong>
</a></li><li class="divider fly-out-top-item"></li>

<li class=""><a title="Incidents" data-qa-selector="operations_incidents_link" href="/gmsh/gmsh/-/incidents"><span>
Incidents
</span>
</a></li><li class=""><a title="Environments" class="shortcuts-environments qa-operations-environments-link" href="/gmsh/gmsh/-/environments"><span>
Environments
</span>
</a></li></ul>
</li>
<li class=""><a data-qa-selector="analytics_anchor" href="/gmsh/gmsh/-/value_stream_analytics"><div class="nav-icon-container">
<svg class="s16" data-testid="chart-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#chart"></use></svg>
</div>
<span class="nav-item-name" data-qa-selector="analytics_link">
Analytics
</span>
</a><ul class="sidebar-sub-level-items" data-qa-selector="analytics_sidebar_submenu">
<li class="fly-out-top-item"><a href="/gmsh/gmsh/-/value_stream_analytics"><strong class="fly-out-top-item-name">
Analytics
</strong>
</a></li><li class="divider fly-out-top-item"></li>
<li class=""><a title="CI / CD" href="/gmsh/gmsh/-/pipelines/charts"><span>CI / CD</span>
</a></li><li class=""><a class="shortcuts-repository-charts" title="Repository" href="/gmsh/gmsh/-/graphs/master/charts"><span>Repository</span>
</a></li><li class=""><a class="shortcuts-project-cycle-analytics" title="Value Stream" href="/gmsh/gmsh/-/value_stream_analytics"><span>Value Stream</span>
</a></li></ul>
</li>
<li class=""><a class="shortcuts-wiki" data-qa-selector="wiki_link" href="/gmsh/gmsh/-/wikis/home"><div class="nav-icon-container">
<svg class="s16" data-testid="book-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#book"></use></svg>
</div>
<span class="nav-item-name">
Wiki
</span>
</a><ul class="sidebar-sub-level-items is-fly-out-only">
<li class="fly-out-top-item"><a href="/gmsh/gmsh/-/wikis/home"><strong class="fly-out-top-item-name">
Wiki
</strong>
</a></li></ul>
</li>
<li class=""><a class="shortcuts-snippets" data-qa-selector="snippets_link" href="/gmsh/gmsh/-/snippets"><div class="nav-icon-container">
<svg class="s16" data-testid="snippet-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#snippet"></use></svg>
</div>
<span class="nav-item-name">
Snippets
</span>
</a><ul class="sidebar-sub-level-items is-fly-out-only">
<li class="fly-out-top-item"><a href="/gmsh/gmsh/-/snippets"><strong class="fly-out-top-item-name">
Snippets
</strong>
</a></li></ul>
</li><li class=""><a title="Members" class="qa-members-link" id="js-onboarding-members-link" href="/gmsh/gmsh/-/project_members"><div class="nav-icon-container">
<svg class="s16" data-testid="users-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#users"></use></svg>
</div>
<span class="nav-item-name">
Members
</span>
</a><ul class="sidebar-sub-level-items is-fly-out-only">
<li class="fly-out-top-item"><a href="/gmsh/gmsh/-/project_members"><strong class="fly-out-top-item-name">
Members
</strong>
</a></li></ul>
</li>
<a class="toggle-sidebar-button js-toggle-sidebar qa-toggle-sidebar rspec-toggle-sidebar" role="button" title="Toggle sidebar" type="button">
<svg class="s16 icon-chevron-double-lg-left" data-testid="chevron-double-lg-left-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#chevron-double-lg-left"></use></svg>
<svg class="s16 icon-chevron-double-lg-right" data-testid="chevron-double-lg-right-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#chevron-double-lg-right"></use></svg>
<span class="collapse-text">Collapse sidebar</span>
</a>
<button name="button" type="button" class="close-nav-button"><svg class="s16" data-testid="close-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#close"></use></svg>
<span class="collapse-text">Close sidebar</span>
</button>
<li class="hidden">
<a title="Activity" class="shortcuts-project-activity" href="/gmsh/gmsh/activity"><span>
Activity
</span>
</a></li>
<li class="hidden">
<a title="Network" class="shortcuts-network" href="/gmsh/gmsh/-/network/master">Graph
</a></li>
<li class="hidden">
<a class="shortcuts-new-issue" href="/gmsh/gmsh/-/issues/new">Create a new issue
</a></li>
<li class="hidden">
<a title="Jobs" class="shortcuts-builds" href="/gmsh/gmsh/-/jobs">Jobs
</a></li>
<li class="hidden">
<a title="Commits" class="shortcuts-commits" href="/gmsh/gmsh/-/commits/master">Commits
</a></li>
<li class="hidden">
<a title="Issue Boards" class="shortcuts-issue-boards" href="/gmsh/gmsh/-/boards">Issue Boards</a>
</li>
</ul>
</div>
</div>

<div class="content-wrapper">
<div class="mobile-overlay"></div>

<div class="alert-wrapper gl-force-block-formatting-context">














<nav class="breadcrumbs container-fluid container-limited" role="navigation">
<div class="breadcrumbs-container">
<button name="button" type="button" class="toggle-mobile-nav"><span class="sr-only">Open sidebar</span>
<svg class="s16" data-testid="hamburger-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#hamburger"></use></svg>
</button><div class="breadcrumbs-links js-title-container" data-qa-selector="breadcrumb_links_content">
<ul class="list-unstyled breadcrumbs-list js-breadcrumbs-list">
<li><a class="group-path breadcrumb-item-text js-breadcrumb-item-text " href="/gmsh"><img class="avatar-tile lazy" width="15" height="15" data-src="/uploads/-/system/group/avatar/2/gmsh-no-text.png" src="data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEKAAEALAAAAAABAAEAAAICTAEAOw==" />gmsh</a><svg class="s8 breadcrumbs-list-angle" data-testid="angle-right-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#angle-right"></use></svg></li> <li><a href="/gmsh/gmsh"><img alt="gmsh" class="avatar-tile lazy" width="15" height="15" data-src="/uploads/-/system/project/avatar/3/gmsh.png" src="data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEKAAEALAAAAAABAAEAAAICTAEAOw==" /><span class="breadcrumb-item-text js-breadcrumb-item-text">gmsh</span></a><svg class="s8 breadcrumbs-list-angle" data-testid="angle-right-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#angle-right"></use></svg></li>

<li>
<h2 class="breadcrumbs-sub-title"><a href="/gmsh/gmsh/blob/master/Common/GmshDefines.h">Repository</a></h2>
</li>
</ul>
</div>
<script type="application/ld+json">
{"@context":"https://schema.org","@type":"BreadcrumbList","itemListElement":[{"@type":"ListItem","position":1,"name":"gmsh","item":"https://gitlab.onelab.info/gmsh"},{"@type":"ListItem","position":2,"name":"gmsh","item":"https://gitlab.onelab.info/gmsh/gmsh"},{"@type":"ListItem","position":3,"name":"Repository","item":"https://gitlab.onelab.info/gmsh/gmsh/blob/master/Common/GmshDefines.h"}]}

</script>

</div>
</nav>

</div>
<div class="container-fluid container-limited ">
<div class="content" id="content-body" itemscope itemtype="http://schema.org/SoftwareSourceCode">
<div class="flash-container flash-container-page sticky" data-qa-selector="flash_container">
</div>

<div class="js-signature-container" data-signatures-path="/gmsh/gmsh/-/commits/29dfc089fb544ca0df3e102f1b0b4db909ccd7dc/signatures?limit=1"></div>

<div class="tree-holder" id="tree-holder">
<div class="nav-block">
<div class="tree-ref-container">
<div class="tree-ref-holder">
<form class="project-refs-form" action="/gmsh/gmsh/-/refs/switch" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" /><input type="hidden" name="destination" id="destination" value="blob" />
<input type="hidden" name="path" id="path" value="Common/GmshDefines.h" />
<div class="dropdown">
<button class="dropdown-menu-toggle js-project-refs-dropdown qa-branches-select" type="button" data-toggle="dropdown" data-selected="master" data-ref="master" data-refs-url="/gmsh/gmsh/refs?sort=updated_desc" data-field-name="ref" data-submit-form-on-click="true" data-visit="true"><span class="dropdown-toggle-text ">master</span><i aria-hidden="true" data-hidden="true" class="fa fa-chevron-down"></i></button>
<div class="dropdown-menu dropdown-menu-paging dropdown-menu-selectable git-revision-dropdown qa-branches-dropdown">
<div class="dropdown-page-one">
<div class="dropdown-title gl-display-flex"><span class="gl-ml-auto">Switch branch/tag</span><button class="dropdown-title-button dropdown-menu-close gl-ml-auto" aria-label="Close" type="button"><svg class="s16 dropdown-menu-close-icon" data-testid="close-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#close"></use></svg></button></div>
<div class="dropdown-input"><input type="search" id="" class="dropdown-input-field qa-dropdown-input-field" placeholder="Search branches and tags" autocomplete="off" /><svg class="s16 dropdown-input-search" data-testid="search-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#search"></use></svg><svg class="s16 dropdown-input-clear js-dropdown-input-clear" data-testid="close-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#close"></use></svg></div>
<div class="dropdown-content"></div>
<div class="dropdown-loading"><div class="gl-spinner-container"><span class="gl-spinner gl-spinner-orange gl-spinner-md gl-mt-7" aria-label="Loading"></span></div></div>
</div>
</div>
</div>
</form>
</div>
<ul class="breadcrumb repo-breadcrumb">
<li class="breadcrumb-item">
<a href="/gmsh/gmsh/-/tree/master">gmsh
</a></li>
<li class="breadcrumb-item">
<a href="/gmsh/gmsh/-/tree/master/Common">Common</a>
</li>
<li class="breadcrumb-item">
<a href="/gmsh/gmsh/-/blob/master/Common/GmshDefines.h"><strong>GmshDefines.h</strong>
</a></li>
</ul>
</div>
<div class="tree-controls gl-children-ml-sm-3"><a class="gl-button btn shortcuts-find-file" rel="nofollow" href="/gmsh/gmsh/-/find_file/master">Find file
</a><a class="gl-button btn js-blob-blame-link" href="/gmsh/gmsh/-/blame/master/Common/GmshDefines.h">Blame</a><a class="gl-button btn" href="/gmsh/gmsh/-/commits/master/Common/GmshDefines.h">History</a><a class="gl-button btn js-data-file-blob-permalink-url" href="/gmsh/gmsh/-/blob/6229c0452e0f4b22450a17ac0e872f32db0c4554/Common/GmshDefines.h">Permalink</a></div>
</div>

<div class="info-well d-none d-sm-block">
<div class="well-segment">
<ul class="blob-commit-info">
<li class="commit flex-row js-toggle-container" id="commit-29dfc089">
<div class="avatar-cell d-none d-sm-block">
<a href="/geuzaine"><img alt="Christophe Geuzaine&#39;s avatar" src="/uploads/-/system/user/avatar/1/avatar.png?width=40" class="avatar s40 d-none d-sm-inline-block" title="Christophe Geuzaine" /></a>
</div>
<div class="commit-detail flex-list">
<div class="commit-content" data-qa-selector="commit_content">
<a class="commit-row-message item-title js-onboarding-commit-item " href="/gmsh/gmsh/-/commit/29dfc089fb544ca0df3e102f1b0b4db909ccd7dc">xmt export</a>
<span class="commit-row-message d-inline d-sm-none">
&middot;
29dfc089
</span>
<div class="committer">
<a class="commit-author-link js-user-link" data-user-id="1" href="/geuzaine">Christophe Geuzaine</a> authored <time class="js-timeago" title="Jun 3, 2020 7:47am" datetime="2020-06-03T07:47:30Z" data-toggle="tooltip" data-placement="bottom" data-container="body">Jun 03, 2020</time>
</div>

</div>
<div class="commit-actions flex-row">

<a class="ci-status-link ci-status-icon-failed d-inline-flex " title="Pipeline: failed" data-toggle="tooltip" data-placement="left" data-container="body" href="/gmsh/gmsh/-/commit/29dfc089fb544ca0df3e102f1b0b4db909ccd7dc/pipelines?ref=master"><svg class="s24" data-testid="status_failed-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#status_failed"></use></svg></a>
<div class="js-commit-pipeline-status" data-endpoint="/gmsh/gmsh/-/commit/29dfc089fb544ca0df3e102f1b0b4db909ccd7dc/pipelines?ref=master"></div>
<div class="commit-sha-group btn-group d-none d-sm-flex">
<div class="label label-monospace monospace">
29dfc089
</div>
<button class="btn gl-button btn btn-default" data-toggle="tooltip" data-placement="bottom" data-container="body" data-title="Copy commit SHA" data-class="gl-button btn btn-default" data-clipboard-text="29dfc089fb544ca0df3e102f1b0b4db909ccd7dc" type="button" title="Copy commit SHA" aria-label="Copy commit SHA"><svg class="s16" data-testid="copy-to-clipboard-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#copy-to-clipboard"></use></svg></button>

</div>
</div>
</div>
</li>

</ul>
</div>


</div>
<div class="blob-content-holder" id="blob-content-holder">
<article class="file-holder">
<div class="js-file-title file-title-flex-parent">
<div class="file-header-content">
<svg class="s16" data-testid="doc-text-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#doc-text"></use></svg>
<strong class="file-title-name gl-word-break-all" data-qa-selector="file_name_content">
GmshDefines.h
</strong>
<button class="btn btn-clipboard btn-transparent" data-toggle="tooltip" data-placement="bottom" data-container="body" data-class="btn-clipboard btn-transparent" data-title="Copy file path" data-clipboard-text="{&quot;text&quot;:&quot;Common/GmshDefines.h&quot;,&quot;gfm&quot;:&quot;`Common/GmshDefines.h`&quot;}" type="button" title="Copy file path" aria-label="Copy file path"><svg class="s16" data-testid="copy-to-clipboard-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#copy-to-clipboard"></use></svg></button>
<small class="mr-1">
6.97 KB
</small>
</div>

<div class="file-actions gl-display-flex gl-flex-fill-1 gl-align-self-start gl-md-justify-content-end"><a class="btn btn-primary js-edit-blob gl-mr-3  btn-sm" data-track-event="click_edit" data-track-label="Edit" href="/gmsh/gmsh/-/edit/master/Common/GmshDefines.h">Edit</a><a class="btn btn-primary ide-edit-button gl-mr-3 btn-inverted btn-sm" data-track-event="click_edit_ide" data-track-label="Web IDE" data-track-property="secondary" href="/-/ide/project/gmsh/gmsh/edit/master/-/Common/GmshDefines.h">Web IDE</a><div class="btn-group ml-2" role="group">

</div><div class="btn-group ml-2" role="group">
<button class="btn btn btn-sm js-copy-blob-source-btn" data-toggle="tooltip" data-placement="bottom" data-container="body" data-class="btn btn-sm js-copy-blob-source-btn" data-title="Copy file contents" data-clipboard-target=".blob-content[data-blob-id=&#39;e93861ed9a61329c40fc572c6155c89aae5cca0f&#39;]" type="button" title="Copy file contents" aria-label="Copy file contents"><svg class="s16" data-testid="copy-to-clipboard-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#copy-to-clipboard"></use></svg></button>
<a class="btn btn-sm has-tooltip" target="_blank" rel="noopener noreferrer" aria-label="Open raw" title="Open raw" data-container="body" href="/gmsh/gmsh/-/raw/master/Common/GmshDefines.h"><svg class="s16" data-testid="doc-code-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#doc-code"></use></svg></a>
<a download="Common/GmshDefines.h" class="btn btn-sm has-tooltip" target="_blank" rel="noopener noreferrer" aria-label="Download" title="Download" data-container="body" href="/gmsh/gmsh/-/raw/master/Common/GmshDefines.h?inline=false"><svg class="s16" data-testid="download-icon"><use xlink:href="/assets/icons-e0a66cb8e6ca64bcdd2a8f111cbd9e94cf727c1bd5939bc71619df5c973fbc87.svg#download"></use></svg></a>

</div></div>
</div>



<div data-blob-data="// Gmsh - Copyright (C) 1997-2020 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef GMSH_DEFINES_H
#define GMSH_DEFINES_H

// IO file formats (numbers should not be changed)
#define FORMAT_MSH          1
#define FORMAT_UNV          2
#define FORMAT_XPM          4
#define FORMAT_PS           5
#define FORMAT_BMP          6
#define FORMAT_GIF          7
#define FORMAT_GEO          8
#define FORMAT_JPEG         9
#define FORMAT_AUTO         10
#define FORMAT_PPM          11
#define FORMAT_YUV          12
#define FORMAT_OPT          15
#define FORMAT_VTK          16
#define FORMAT_MPEG         17
#define FORMAT_TEX          18
#define FORMAT_VRML         19
#define FORMAT_EPS          20
#define FORMAT_MAIL         21
#define FORMAT_PNG          22
#define FORMAT_TXT          23
#define FORMAT_PDF          24
#define FORMAT_RMED         25
#define FORMAT_POS          26
#define FORMAT_STL          27
#define FORMAT_P3D          28
#define FORMAT_SVG          29
#define FORMAT_MESH         30
#define FORMAT_BDF          31
#define FORMAT_CGNS         32
#define FORMAT_MED          33
#define FORMAT_DIFF         34
#define FORMAT_BREP         35
#define FORMAT_STEP         36
#define FORMAT_IGES         37
#define FORMAT_IR3          38
#define FORMAT_INP          39
#define FORMAT_PLY2         40
#define FORMAT_CELUM        41
#define FORMAT_SU2          42
#define FORMAT_MPEG_PREVIEW 43
#define FORMAT_PGF          44
#define FORMAT_PVTU         45
#define FORMAT_X3D          46
#define FORMAT_TOCHNOG      47
#define FORMAT_TIKZ         48
#define FORMAT_NEU          49
#define FORMAT_MATLAB       50
#define FORMAT_KEY          51
#define FORMAT_XMT          52

// Element types
#define TYPE_PNT     1
#define TYPE_LIN     2
#define TYPE_TRI     3
#define TYPE_QUA     4
#define TYPE_TET     5
#define TYPE_PYR     6
#define TYPE_PRI     7
#define TYPE_HEX     8
#define TYPE_POLYG   9
#define TYPE_POLYH   10
#define TYPE_XFEM    11
#define TYPE_MINI    12
#define TYPE_TRIH    13
#define TYPE_MAX_NUM 13 // keep this up-to-date when adding new type

// Element types in .msh file format (numbers should not be changed)
#define MSH_LIN_2    1
#define MSH_TRI_3    2
#define MSH_QUA_4    3
#define MSH_TET_4    4
#define MSH_HEX_8    5
#define MSH_PRI_6    6
#define MSH_PYR_5    7
#define MSH_LIN_3    8
#define MSH_TRI_6    9
#define MSH_QUA_9    10
#define MSH_TET_10   11
#define MSH_HEX_27   12
#define MSH_PRI_18   13
#define MSH_PYR_14   14
#define MSH_PNT      15
#define MSH_QUA_8    16
#define MSH_HEX_20   17
#define MSH_PRI_15   18
#define MSH_PYR_13   19
#define MSH_TRI_9    20
#define MSH_TRI_10   21
#define MSH_TRI_12   22
#define MSH_TRI_15   23
#define MSH_TRI_15I  24
#define MSH_TRI_21   25
#define MSH_LIN_4    26
#define MSH_LIN_5    27
#define MSH_LIN_6    28
#define MSH_TET_20   29
#define MSH_TET_35   30
#define MSH_TET_56   31
#define MSH_TET_22   32
#define MSH_TET_28   33
#define MSH_POLYG_   34
#define MSH_POLYH_   35
#define MSH_QUA_16   36
#define MSH_QUA_25   37
#define MSH_QUA_36   38
#define MSH_QUA_12   39
#define MSH_QUA_16I  40
#define MSH_QUA_20   41
#define MSH_TRI_28   42
#define MSH_TRI_36   43
#define MSH_TRI_45   44
#define MSH_TRI_55   45
#define MSH_TRI_66   46
#define MSH_QUA_49   47
#define MSH_QUA_64   48
#define MSH_QUA_81   49
#define MSH_QUA_100  50
#define MSH_QUA_121  51
#define MSH_TRI_18   52
#define MSH_TRI_21I  53
#define MSH_TRI_24   54
#define MSH_TRI_27   55
#define MSH_TRI_30   56
#define MSH_QUA_24   57
#define MSH_QUA_28   58
#define MSH_QUA_32   59
#define MSH_QUA_36I  60
#define MSH_QUA_40   61
#define MSH_LIN_7    62
#define MSH_LIN_8    63
#define MSH_LIN_9    64
#define MSH_LIN_10   65
#define MSH_LIN_11   66
#define MSH_LIN_B    67
#define MSH_TRI_B    68
#define MSH_POLYG_B  69
#define MSH_LIN_C    70
// TETS COMPLETE (6-&gt;10)
#define MSH_TET_84   71
#define MSH_TET_120  72
#define MSH_TET_165  73
#define MSH_TET_220  74
#define MSH_TET_286  75
// TETS INCOMPLETE (6-&gt;10)
#define MSH_TET_34   79
#define MSH_TET_40   80
#define MSH_TET_46   81
#define MSH_TET_52   82
#define MSH_TET_58   83
//
#define MSH_LIN_1    84
#define MSH_TRI_1    85
#define MSH_QUA_1    86
#define MSH_TET_1    87
#define MSH_HEX_1    88
#define MSH_PRI_1    89
#define MSH_PRI_40   90
#define MSH_PRI_75   91
// HEXES COMPLETE (3-&gt;9)
#define MSH_HEX_64   92
#define MSH_HEX_125  93
#define MSH_HEX_216  94
#define MSH_HEX_343  95
#define MSH_HEX_512  96
#define MSH_HEX_729  97
#define MSH_HEX_1000 98
// HEXES INCOMPLETE (3-&gt;9)
#define MSH_HEX_32   99
#define MSH_HEX_44   100
#define MSH_HEX_56   101
#define MSH_HEX_68   102
#define MSH_HEX_80   103
#define MSH_HEX_92   104
#define MSH_HEX_104  105
// PRISMS COMPLETE (5-&gt;9)
#define MSH_PRI_126  106
#define MSH_PRI_196  107
#define MSH_PRI_288  108
#define MSH_PRI_405  109
#define MSH_PRI_550  110
// PRISMS INCOMPLETE (3-&gt;9)
#define MSH_PRI_24   111
#define MSH_PRI_33   112
#define MSH_PRI_42   113
#define MSH_PRI_51   114
#define MSH_PRI_60   115
#define MSH_PRI_69   116
#define MSH_PRI_78   117
// PYRAMIDS COMPLETE (3-&gt;9)
#define MSH_PYR_30   118
#define MSH_PYR_55   119
#define MSH_PYR_91   120
#define MSH_PYR_140  121
#define MSH_PYR_204  122
#define MSH_PYR_285  123
#define MSH_PYR_385  124
// PYRAMIDS INCOMPLETE (3-&gt;9)
#define MSH_PYR_21   125
#define MSH_PYR_29   126
#define MSH_PYR_37   127
#define MSH_PYR_45   128
#define MSH_PYR_53   129
#define MSH_PYR_61   130
#define MSH_PYR_69   131
// Additional types
#define MSH_PYR_1    132
#define MSH_PNT_SUB  133
#define MSH_LIN_SUB  134
#define MSH_TRI_SUB  135
#define MSH_TET_SUB  136
#define MSH_TET_16   137
#define MSH_TRI_MINI 138
#define MSH_TET_MINI 139
#define MSH_TRIH_4   140
#define MSH_MAX_NUM  140 // keep this up-to-date when adding new type

// Geometric entities
#define ENT_NONE    0
#define ENT_POINT   (1&lt;&lt;0)
#define ENT_CURVE   (1&lt;&lt;1)
#define ENT_SURFACE (1&lt;&lt;2)
#define ENT_VOLUME  (1&lt;&lt;3)
#define ENT_ALL     (ENT_POINT | ENT_CURVE | ENT_SURFACE | ENT_VOLUME)

// 2D meshing algorithms (numbers should not be changed)
#define ALGO_2D_MESHADAPT         1
#define ALGO_2D_AUTO              2
#define ALGO_2D_INITIAL_ONLY      3
#define ALGO_2D_DELAUNAY          5
#define ALGO_2D_FRONTAL           6
#define ALGO_2D_BAMG              7
#define ALGO_2D_FRONTAL_QUAD      8
#define ALGO_2D_PACK_PRLGRMS      9
#define ALGO_2D_PACK_PRLGRMS_CSTR 10

// 3D meshing algorithms (numbers should not be changed)
#define ALGO_3D_DELAUNAY     1
#define ALGO_3D_INITIAL_ONLY 3
#define ALGO_3D_FRONTAL      4
#define ALGO_3D_MMG3D        7
#define ALGO_3D_RTREE        9
#define ALGO_3D_HXT          10

// Meshing methods
#define MESH_NONE         0
#define MESH_TRANSFINITE  1
#define MESH_UNSTRUCTURED 2

// QuadTri options (structured/unstructured coupling with pyramids)
#define NO_QUADTRI                0
#define QUADTRI_ADDVERTS_1        1
#define QUADTRI_ADDVERTS_1_RECOMB 2
#define QUADTRI_NOVERTS_1         3
#define QUADTRI_NOVERTS_1_RECOMB  4
#define TRANSFINITE_QUADTRI_1     5

#endif
" data-is-ci-config-file="false" id="js-blob-toggle-graph-preview"></div>
<div class="blob-viewer" data-path="Common/GmshDefines.h" data-type="simple" data-url="/gmsh/gmsh/-/blob/master/Common/GmshDefines.h?format=json&amp;viewer=simple">
<div class="text-center gl-mt-4 gl-mb-3">
<span class="gl-spinner gl-spinner-orange gl-spinner-md qa-spinner" aria-label="Loading"></span>
</div>

</div>


</article>
</div>

<div class="modal" id="modal-upload-blob">
<div class="modal-dialog modal-lg">
<div class="modal-content">
<div class="modal-header">
<h3 class="page-title">Replace GmshDefines.h</h3>
<button aria-label="Close" class="close" data-dismiss="modal" type="button">
<span aria-hidden>&times;</span>
</button>
</div>
<div class="modal-body">
<form class="js-quick-submit js-upload-blob-form" data-method="put" action="/gmsh/gmsh/-/update/master/Common/GmshDefines.h" accept-charset="UTF-8" method="post"><input name="utf8" type="hidden" value="&#x2713;" /><input type="hidden" name="_method" value="put" /><input type="hidden" name="authenticity_token" value="noywC6+dS7e0vSBqokI9PYXbXI1ietXO0T7WflXHBU/+RvGcONZAdvQwtIfcw2OENx1NmCjwe9R60AX0+1G0ow==" /><div class="dropzone">
<div class="dropzone-previews blob-upload-dropzone-previews">
<p class="dz-message light">
Attach a file by drag &amp; drop or <a class="markdown-selector" href="#">click to upload</a>
</p>
</div>
</div>
<br>
<div class="dropzone-alerts gl-alert gl-alert-danger gl-mb-5 data" style="display:none"></div>
<div class="form-group row commit_message-group">
<label class="col-form-label col-sm-2" for="commit_message-6a88fb4dbfff582c55a2de21e2f70095">Commit message
</label><div class="col-sm-10">
<div class="commit-message-container">
<div class="max-width-marker"></div>
<textarea name="commit_message" id="commit_message-6a88fb4dbfff582c55a2de21e2f70095" class="form-control js-commit-message" placeholder="Replace GmshDefines.h" required="required" rows="3">
Replace GmshDefines.h</textarea>
</div>
</div>
</div>

<input type="hidden" name="branch_name" id="branch_name" />
<input type="hidden" name="create_merge_request" id="create_merge_request" value="1" />
<input type="hidden" name="original_branch" id="original_branch" value="master" class="js-original-branch" />

<div class="form-actions">
<button name="button" type="button" class="btn gl-button btn-success btn-upload-file" id="submit-all"><div class="spinner spinner-sm gl-mr-2 js-loading-icon hidden"></div>
Replace file
</button><a class="btn gl-button btn-cancel" data-dismiss="modal" href="#">Cancel</a>
<div class="inline gl-ml-3">
A new branch will be created in your fork and a new merge request will be started.
</div>

</div>
</form></div>
</div>
</div>
</div>

</div>


</div>
</div>
</div>
</div>


<script>
//<![CDATA[
if ('loading' in HTMLImageElement.prototype) {
  document.querySelectorAll('img.lazy').forEach(img => {
    img.loading = 'lazy';
    let imgUrl = img.dataset.src;
    // Only adding width + height for avatars for now
    if (imgUrl.indexOf('/avatar/') > -1 && imgUrl.indexOf('?') === -1) {
      const targetWidth = img.getAttribute('width') || img.width;
      imgUrl += `?width=${targetWidth}`;
    }
    img.src = imgUrl;
    img.removeAttribute('data-src');
    img.classList.remove('lazy');
    img.classList.add('js-lazy-loaded', 'qa-js-lazy-loaded');
  });
}

//]]>
</script>

</body>
</html>

