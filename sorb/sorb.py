import os
import sys

#("shakespeare","","")
#("bible","","")

sorb_extracts = [("shakespeare", "A beggar's book out-worths a noble's blood.", "King Henry VIII, 1:1"),
("shakespeare", "When I had lost one shaft, I shot his fellow of the selfsame flight, the selfsame way, with more advised watch.", "Merchant of Venice 1:1"),
("bible", "Who hath measured the waters in the hollow of his hand and meted out heaven with the span, and comprehended the dust of the earth in a measure, and weighted the mountains in scales, and the hills in a balance?", "Isaiah 40:12"),
("wodehouse", "It was my Uncle George who discovered that alcohol was a food well in advance of modern medical thought.","The Inimitable Jeeves"),
("shakespeare", "What a devil hast thou to do with the time of the day?", "King Henry The Fourth part i 1:1"),
("shakespeare", "Thou hast called her to a reckoning many a time and oft.", "King Henry The Fourth part i 1:1"),
("shakespeare", "Wisdom cries out in the streets and no man regards it.", "King Henry The Fourth part i 1:1"),
("bible", "Wisdom crieth without; she uttereth her voice in the streets", "Proverbs 1:20"),
("bible", "Go to the ant, thou sluggard ; consider her ways, and be wise.", "Proverbs 6:6"),
("bible", "I wisdom dwell with prudence, and find out knowledge of witty inventions.", "Proverbs 8:12"),
("shakespeare", "Disparage not the faith thou dost not know, Lest to thy peril thou aby it...", "A Midsummer Night's Dream 3:2"),
("shakespeare", "Alack the heavy day, That I have worn so many winters out, And know not now what name to call myself.","King Richard The Second 4:1"),
("shakespeare", "For what can we bequeath save our bodies to the ground", "King Richard The Second 3:2"),
("shakespeare", "Cover your heads, and mock not flesh and blood", "King Richard The Second 3:2"),
("shakespeare", "He is a very serpent in my way, And wheresoe'er this foot of mine doth tread He lies before me", "King John 3:3"),
("shakespeare", "Make dust our paper, and with rainy eyes write sorrow on the bosom of the earth.", "King Richard II 3:3"),
("shakespeare", "Bring thou her husband: This is the hole where Aaron bid us hide him.","Titus Andronichus 2:3"),
("shakespeare", "Therefore I say again, I utterly abhor, yea, from my soul Refuse you for my judge.","King Henry VIII 2:4"),
("shakespeare", "There is a vice that most I do abhor, And most desire should meet the blow of justice.","Measure for Measure"),
("shakespeare", "A true devoted pilgrim is not weary To measure kingdoms with his feeble steps.","Three Gentleman of Verona 2:7"),
("shakespeare", "I took by th'throat the circumcised dog, And smote him", "Othello 5:2"),
("shakespeare", "Does your worship mean to geld and splay all the youth of the city?", "Measure for Measure 2:1"),
("shakespeare","`Thou shalt not steal'?", "Measure For Measure 1:2"),
("bible", "Remove not the ancient landmark, which thy fathers have set.", "Proverbs 22:28"),
("shakespeare", "Graves at my command Have waked their sleepers", "Tempest 5:1"),
("shakespeare", "Sleepest or wakest thou jolly shepherd? Thy sheep be in the corn", "King Lear 3:6"),
("shakespeare", "What is a man, If his chief good and market of his time Be but to sleep and feed?","Hamlet 4:4"),
("shakespeare", "The fairest votary took up that fire Which many legions of true hearts had warmed.", "Sonnet 154"),
("shakespeare", "My thoughts, from far where I abide, intend a zealous pilgrimage to thee.", "Sonnet 27"),
("bible", "She openeth her mouth with wisdom; and in her tongue is the law of kindness.", "Proverbs 31:26"),
("bible","Foolishness is bound in the heart of a child; but the rod of correction shall drive it far from him.", "Proverbs 22:15"),
("bible", "Pride goeth before destruction, and a haughty spirit before a fall.", "Proverbs 16:18"),
("bible", "The glory of young men is their strength: and the beauty of old men is their gray head.","Proverbs 20:26"),
("bible", "A man hath joy by the answer of his mouth: and a word spoken in due season, how good is it!", "Proverbs 15:23"),
("bible", "Yet a little sleep, a little slumber, a little folding of the hands to sleep","Proverbs 6:10"),
("bible","Thou shalt not steal.", "Exodus 20:15"),
("bible","A fugitive and a vagabond shalt thou be in the earth.","Genesis 4:12"),
("bible","Wherefore is this noise of the city being in an uproar?", "Kings I 1:41"),
("bible","I desired Titus, and with him I sent a brother. Did Titus make a gain of you? Walked we not in the same spirit? Walked we not in the same steps?","Corinthians II 12:18"),
("bible","How doth the city sit solitary, that was full of people! How is she become as a widow!","Lamentations 1:1"),
("bible", "All my inward friends abhorred me: and they whom I loved are turned against me.", "Job 19:18"),
("bible", "It is like the precious ointment upon the head, that ran down upon the beard, even Aaron's beard","Psalms 133:2"),
("bible", "Every man at the beginning doth set forth good wine; and when men have well drunk, then has which is worse.", "John 2:10"),
("bible", "It is better to dwell in the wilderness, than with a contentious and an angry woman.", "Proverbs 21:19"),
("bible", "Give strong drink unto him that is ready to perish, and wine unto those that be of heavy hearts.", "Proverbs 31:6"),
("bible", "Let him drink, and forget his poverty, and remember his misery no more.", "Proverbs 31:7"),
("bible", "She is not afraid of the snow for her household: for all her household are clothed with scarlet.", "Proverbs 31:21"),
("bible", "they are like the deaf adder that stoppeth her ear; Which will not hearken to the voice of charmers, charming never so wisely.", "Psalm 58:4"),
("bible", "Lo, this is the man that made not God his strength; but trusted in the abundance of his riches, and strengthened himself in his wickedness","Psalm 52:7"),
("bible", "thou art beside thyself; much learning doth make thee mad", "Acts 26:24-25"),
("bible", "Surely in vain the net is spread in the sight of any bird.", "Proverbs 1:17"),
("bible", "A man after his own heart", "Samuel 13:14"),
("bible", "There is a path which no foul knoweth and which the vulture's eye has not seen.", "Job 28:7"),
("bible", "Thou compassest my path and my lying down, and art acquainted with all my ways.","Psalms 139:3"),
("bible", "If I take the wings of the morning, and dwell in the uttermost parts of the sea...","Psalms 139:9", "... Even there shall thy hand lead me, and thy right hand shall hold me."),
("bible", "I was made in secret, and curiously wrought in the lowest part of the earth.", "Psalms 139:15"),
("bible", "There is no new thing under the sun.","Ecclesiastes 1:9"),
("bible", "For in much wisdom is much grief: and he that increaseth knowledge increaseth sorrow.", "Ecclesiastes 18:1"),
("bible", "Who can find a virtuous woman? for her price is far above rubies.","Proverbs 31:10"),
("bible", "I go the way of all the earth, be thou strong therefore, and show thyself a man.","I Kings 2:2"),
("bible", " I am made all things to all men, that I might by all means save some.", "I Corinthians 9:22"),
("wodehouse", "... the children who made mock of the prophet Elisha.","Introduction to Summer Lightning"),
("wodehouse", "Unseen, in the background, Fate was quietly slipping the lead into the boxing-glove.","Very Good, Jeeves"),
("wodehouse", "We Scripture-knowledge sharks stick together.", "Right Ho, Jeeves"),
("wodehouse", "It was impossible for him to have won the Scripture-knowledge prize without systematic cheating on an impressive scale.",
"Right Ho, Jeeves"),
("bible", "I have been a stranger in a strange land.", "Exodus 2:22"),
("bible", "Because there were no graves in Egypt, hast thou taken us away to die in the wilderness?", "Exodus 14:11"),
("bible", "he hath triumphed gloriously: the horse and his rider hath he thrown into the sea.","Exodus 15:1"),
("shakespeare", "He hath the horn of abundance, and the lightness of his wife shines through it: and yet cannot he see, though he have his own lantern to him", "Henry iv part 2 1:2")
]
