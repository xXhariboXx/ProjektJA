   M S B A M L   `   `   `  ����J  FWpfApplication1, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null1clr-namespace:WpfApplication1WpfApplication1  5   	        WpfApplication1.MainWindow   5   	   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    5   	   R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     5   	   @d2http://schemas.microsoft.com/expression/blend/2008     5   	   Jmc;http://schemas.openxmlformats.org/markup-compatibility/2006     5   	   /localclr-namespace:WpfApplication1     5   	   -   5        $��Main_window��5   	     � Title$  
MainWindow��5	   	   $��437.967��6   $��716.443��6-   .��5
      � $
��Left=�6   $	��400��6&   $��-4,1,0,0q�63   $	��Top=�6E   $	��686��6]    � ColumnDefinitions 5   
   �� 5      $
��178* �6   �� 5      $
��165* �6   .+�5   
   �� -   5         ��Name$
 Btn1��5      $��
Zatwierdź��6   $
��Left=�63   $��23��6N   $��
177,73,0,0q�6Z   $	��Top=�6n   $��65��6�   6�   �� 5   
   -   5         ��Name$ txtBox1��5      $
��Left=�6#   $��24��6>   $��	37,72,0,0q�6J    �� TextWrapping$
 Wrap=�6]   $	��Top=�6q   $	��118��6�   �� 5   
   -   5         ��Name$ dllASM_chbx��5      $��DLL ASM��6(   $
��Left=�6:   $��19��6U   $��
37,227,0,0q�6a   $	��Top=�6u   $	��118��6�   �� 5   
   -   5        $ 	dllC_chbx��5      $��DLL C#��6&   $
��Left=�67   $��19��6R   $��
37,246,0,0q�6^   $	��Top=�6r   $	��118��6�   �� 5   
   -   5        $ dllCOMP_chbx��5      $��Test porównawczy��6)   $
��Left=�6D   $��19��6_   $��
37,265,0,0q�6k   $	��Top=�6   $	��118��6�   &   System.Windows.Controls.DataGrid  5   
   -   5          Name$ equations_datagrid��5      $
��Left=�6/   $��
316,42,0,0q�6J   $	��Top=�6^   $
��Auto��6v   $
��Auto��6�   $��2��6�   $��0q�6�      CanUserAddRows .     CanUserDeleteRows .  	   SelectionUnit$
	 Cell=�6�   
   EnableRowVirtualization
 .   J� IsManipulationEnabled .    AutoGenerateColumns .  $y�Hidden=�6X  �� 5   
   -   5        $ calculate_btn��5      $
��Licz��6(   $
��Left=�67   $��23��6R   $��223.053,327,0,0q�6^   $	��Top=�6w   $��65��6�   6�   �� 5   
   -	   5         ��Name$
 lbl1��5      $ ��Podaj liczbę niewiadomych��6   $
��Left=�6B   $��25��6]   $��	37,42,0,0q�6i   $	��Top=�6|   $	��205��6�     5   
   -
   5        $ result_datagrid��5      $
��Left=�6,   $��316,315,0,0q�6G   $	��Top=�6\   $
��Auto��6t   $
��Auto��6�   $��2��6�   $��0q�6�    .   .  $
	 Cell=�6�   
 .   .  .  $y�Hidden=�6V  �� 5   
   -   5        $ result_label��5      $��Wynik��6&   $��1��66   $
��Left=�6F   $��	0,282,0,0q�6a   $	��Top=�6t   $��71��6�   $y�Hidden=�6�   5      5      