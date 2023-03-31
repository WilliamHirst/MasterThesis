# $Id: russian.perl,v 1.1 2001/04/18 01:39:21 RRM Exp $
#
# russian.perl for russian babel, inspired heavily by german.perl
# by Georgy Salnikov <sge@nmr.nioch.nsc.ru>


package russian;

print " [russian]";

sub main'russian_translation { @_[0] }


package main;

if (defined &addto_languages) { &addto_languages('russian') };

sub russian_titles {
    $toc_title = "����������";
    $lof_title = "������ �����������";
    $lot_title = "������ ������";
    $idx_title = "���������� ���������";
    $ref_title = "������ ����������";
    $bib_title = "����������";
    $abs_title = "���������";
    $app_title = "����������";
    $pre_title = "�����������";
    $foot_title = "����������";
    $thm_title = "�������";
    $fig_name = "���.";
    $tab_name = "�������";
    $prf_name = "��������������";
    $date_name = "����";
    $page_name = "�.";
  #  Sectioning-level titles
    $part_name = "�����";
    $chapter_name = "�����";
    $section_name = "������";
    $subsection_name = "���������";
    $subsubsection_name = "������������";
    $paragraph_name = "��������";
  #  Misc. strings
    $child_name = "���������";
    $info_title = "�� ���� ��������� ...";
    $also_name = "��. �����";
    $see_name = "��.";
  #  names in navigation panels
    $next_name = "����.";
    $up_name = "����";
    $prev_name = "����.";
  #  field names in email headers
    $encl_name = "���.";
    $headto_name = "��.";
    $cc_name = "���.";

    @Month = ('', '������', '�������', '�����', '������', '���',
	      '����', '����', '�������', '��������', '�������',
	      '������', '�������');
    $GENERIC_WORDS = "and|the|of|for|by|a|an|to";
}


sub russian_today {
    local($today) = &get_date();
    $today =~ s|(\d+)/0?(\d+)/|$2 $Month[$1] |;
    join('',$today,$_[0],' �.');
}



# use'em
&russian_titles;
$default_language = 'russian';
$TITLES_LANGUAGE = 'russian';
$russian_encoding = 'koi8-r';

# $Log: russian.perl,v $
# Revision 1.1  2001/04/18 01:39:21  RRM
#      support for the Russian language using KOI8-R encoding, and as an
#      option to the Babel package.
#      supplied by:  Georgy Salnikov  <sge@nmr.nioch.nsc.ru>
#
# Revision 1.1  1998/08/25 02:11:25  RRM
# 	Babel language support
#
#

1;
