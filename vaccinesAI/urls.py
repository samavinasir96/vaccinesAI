"""
URL configuration for vaccinesAI project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path

from bvacai import views as bvacai_views
from eskapeml import views as eskapeml_views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('',eskapeml_views.home,name='index'),

    # path('bvacai_upload/',bvacai_views.upload_sequence,name='bvacai_upload'),
    # path('bvacai_results/',bvacai_views.get_results,name='bvacai_results'),
    # path('bvacai_annotate/', bvacai_views.calculate_features, name='bvacai_annotate_features'),
    # path('bvacai_get_latest_logs/', bvacai_views.get_latest_logs, name='get_latest_logs'),
    # path('bvacai_download_csv/', bvacai_views.download_csv, name='bvacai_download_csv'),
    # path('bvacai_delete_files/', bvacai_views.delete_files, name='bvacai_delete_files'),
    # path('bvacai_display/', bvacai_views.display_features, name='bvacai_display_features'),

    path('get_progress/', eskapeml_views.get_progress, name='get_progress'),

    path('eskapeml_upload/',eskapeml_views.upload_sequence,name='eskapeml_upload'),
    path('eskapeml_results/',eskapeml_views.get_results,name='eskapeml_results'),
    path('eskapeml_annotate/', eskapeml_views.calculate_features, name='eskapeml_annotate_features'),
    path('eskapeml_get_latest_logs/', eskapeml_views.get_latest_logs, name='get_latest_logs'),
    path('eskapeml_download_csv/', eskapeml_views.download_csv, name='eskapeml_download_csv'),
    path('eskapeml_delete_files/', eskapeml_views.delete_files, name='eskapeml_delete_files'),
    path('eskapeml_display/', eskapeml_views.display_features, name='eskapeml_display_features'),
    # path('processing/', eskapeml_views.processing_view, name='processing'),


    path('glossary/',eskapeml_views.glossary,name='glossary'),
    path('faqs/',eskapeml_views.faqs,name='faqs')


]
