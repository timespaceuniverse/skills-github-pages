---
layout: default
title: Blogs on "aha!" moment
---

# Blogs on "Aha!" moments

## Recent Posts
<ul>
  {% for post in site.posts %}
    <li>
      <span style="color: #666; font-family: monospace;">{{ post.date | date: "%Y-%m-%d" }}</span> — 
      <a href="{{ post.url | relative_url }}">{{ post.title }}</a>
    </li>
  {% endfor %}
</ul>
