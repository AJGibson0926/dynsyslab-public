#include "dynsyslab_compute/Workbench.hpp"

#include "dynsyslab_compute/PanelCommon.hpp"

#include <algorithm>

namespace dsl::compute {

void Workbench::layout() {
    RECT area{}; GetClientRect(window_,&area); const int w=area.right,h=area.bottom;
    const auto ui=dsl::theme::metrics(window_);
    const int margin=ui.margin,gap=ui.gap;
    const int header=ui.headerHeight,footer=scalePixels(window_,34);
    const int splitterWidth=scalePixels(window_,8);
    if(!paneWidthsInitialized_){
        navigationWidth_=scalePixels(window_,260);
        inspectorWidth_=scalePixels(window_,300);
        paneWidthsInitialized_=true;
    }
    dsl::HorizontalPaneConstraints constraints;
    constraints.totalWidth=w;constraints.margin=margin;constraints.gap=gap;
    constraints.desiredNavigationWidth=navigationWidth_;
    constraints.desiredInspectorWidth=inspectorWidth_;
    constraints.minimumNavigationWidth=scalePixels(window_,190);
    constraints.maximumNavigationWidth=scalePixels(window_,480);
    constraints.minimumCenterWidth=scalePixels(window_,650);
    constraints.minimumInspectorWidth=scalePixels(window_,240);
    constraints.maximumInspectorWidth=scalePixels(window_,480);
    constraints.inspectorThreshold=scalePixels(window_,1120);
    const auto panes=dsl::resolveHorizontalPaneLayout(constraints);
    const int navW=panes.navigationWidth,inspectorW=panes.inspectorWidth;
    const int centerX=panes.centerX,centerW=panes.centerWidth;
    const int navigationToolsY=header;
    const int toolY=scalePixels(window_,8),toolbarH=ui.toolbarHeight;
    const int refreshW=scalePixels(window_,32);
    if(!solverActive_){
        MoveWindow(search_,margin,navigationToolsY,
            navW-refreshW-scalePixels(window_,6),toolbarH,TRUE);
        MoveWindow(refresh_,margin+navW-refreshW,navigationToolsY,
            refreshW,toolbarH,TRUE);
    }
    const int toolGap=scalePixels(window_,6);
    int toolRight=w-margin;
    auto placeTool=[&](HWND control,int logicalWidth){
        const int width=scalePixels(window_,logicalWidth);
        toolRight-=width;MoveWindow(control,toolRight,toolY,width,toolbarH,TRUE);
        toolRight-=toolGap;
    };
    placeTool(help_,40);placeTool(importData_,96);placeTool(archiveTools_,102);
    placeTool(queue_,96);toolRight-=scalePixels(window_,10);
    placeTool(openJob_,88);placeTool(saveJob_,88);
    if(solverActive_)solveJobs_.layout(margin,header,navW,
        h-footer-margin-header);
    else {
        const int navigationY=navigationToolsY+toolbarH+scalePixels(window_,6);
        MoveWindow(navigation_.handle(),margin,navigationY,navW,
            h-navigationY-footer-margin,TRUE);
    }
    leftSplitter_.move(panes.leftSplitterCenter-splitterWidth/2,header,splitterWidth,
        h-header-footer-margin);
    const int tabH=scalePixels(window_,34),panelGap=scalePixels(window_,8);
    const int solverTabW=scalePixels(window_,118),postTabW=scalePixels(window_,166);
    MoveWindow(solverTab_,centerX,header,solverTabW,tabH,TRUE);
    MoveWindow(postTab_,centerX+solverTabW+scalePixels(window_,7),header,postTabW,tabH,TRUE);
    const int viewerW=scalePixels(window_,124);
    MoveWindow(viewer_,centerX+centerW-viewerW,header,viewerW,tabH,TRUE);
    const int workflowY=header+tabH+panelGap;
    const int workflowH=toolbarH;
    HWND workflowStages[]={stageSystem_,stageInitialConditions_,stageControl_,
        stageIntegrator_,stageReview_};
    const int workflowGap=scalePixels(window_,6);
    const int workflowW=std::max(1,(centerW-workflowGap*4)/5);
    for(int index=0;index<5;++index){
        MoveWindow(workflowStages[index],centerX+index*(workflowW+workflowGap),
            workflowY,workflowW,workflowH,TRUE);
        ShowWindow(workflowStages[index],solverActive_?SW_SHOW:SW_HIDE);
    }
    const int panelY=solverActive_
        ?workflowY+workflowH+panelGap:header+tabH+panelGap;
    const int actionCaptionH=scalePixels(window_,15),actionH=scalePixels(window_,38);
    const int actionGap=scalePixels(window_,7);
    const int actionY=h-footer-actionH;
    const int collapsedReviewH=panes.inspectorVisible?0:scalePixels(window_,160);
    const int collapsedReviewGap=collapsedReviewH?scalePixels(window_,8):0;
    const int contextHeight=std::max(scalePixels(window_,180),
        actionY-actionGap-actionCaptionH-panelY-collapsedReviewH-collapsedReviewGap);
    const int executionGap=scalePixels(window_,10);
    const int executionH=std::min(scalePixels(window_,96),
        std::max(scalePixels(window_,84),contextHeight*22/100));
    const int editorH=std::max(scalePixels(window_,80),
        contextHeight-executionGap-executionH);
    solver_.layoutEditor(centerX,panelY,centerW,editorH);
    solveContext_.layout(centerX,panelY,centerW,editorH);
    solver_.layoutExecution(centerX,panelY+editorH+executionGap,
        centerW,executionH);
    MoveWindow(jobActionsCaption_,centerX,actionY-actionGap-actionCaptionH,
        centerW,actionCaptionH,TRUE);
    solver_.layoutJobActions(centerX,actionY,centerW,actionH);
    if(postCreated_)post_.layout(centerX,panelY,centerW,
        h-panelY-footer-collapsedReviewH-collapsedReviewGap);
    rightSplitter_.show(panes.inspectorVisible);
    if(panes.inspectorVisible){
        ShowWindow(reviewRailTitle_,SW_SHOW);ShowWindow(inspector_.handle(),SW_SHOW);
        ShowWindow(output_,SW_SHOW);const int right=panes.inspectorX;
        rightSplitter_.move(panes.rightSplitterCenter-splitterWidth/2,header,
            splitterWidth,h-header-footer-margin);
        const int rightHeight=h-header-footer;
        const int reviewTitleH=scalePixels(window_,24),reviewGap=scalePixels(window_,8);
        MoveWindow(reviewRailTitle_,right,header,inspectorW,reviewTitleH,TRUE);
        const int railContentY=header+reviewTitleH+reviewGap;
        const int railContentH=rightHeight-reviewTitleH-reviewGap;
        const int inspectorH=railContentH*36/100;
        MoveWindow(inspector_.handle(),right,railContentY,inspectorW,inspectorH,TRUE);
        MoveWindow(output_,right,railContentY+inspectorH+reviewGap,inspectorW,
            railContentH-inspectorH-reviewGap,TRUE);
    } else { ShowWindow(reviewRailTitle_,SW_HIDE);ShowWindow(inspector_.handle(),SW_HIDE);
        ShowWindow(output_,SW_SHOW);
        const int outputY=panelY+contextHeight+collapsedReviewGap;
        MoveWindow(output_,centerX,outputY,centerW,collapsedReviewH,TRUE); }
    const int progressW=scalePixels(window_,220),labelW=scalePixels(window_,105);
    const int cancelW=scalePixels(window_,112),footerPad=scalePixels(window_,7);
    MoveWindow(status_,centerX,h-footer+footerPad,
        std::max(scalePixels(window_,80),w-centerX-margin-progressW-labelW-cancelW-scalePixels(window_,21)),
        footer-scalePixels(window_,9),TRUE);
    MoveWindow(cancelTask_,w-margin-progressW-labelW-cancelW-scalePixels(window_,14),
        h-footer+scalePixels(window_,3),cancelW,scalePixels(window_,26),TRUE);
    MoveWindow(progressLabel_,w-margin-progressW-labelW-scalePixels(window_,7),
        h-footer+footerPad,labelW,footer-scalePixels(window_,9),TRUE);
    MoveWindow(progress_,w-margin-progressW,h-footer+scalePixels(window_,8),
        progressW,scalePixels(window_,18),TRUE);
    InvalidateRect(window_,nullptr,FALSE);
}

void Workbench::paint(HDC dc) {
    RECT area{}; GetClientRect(window_,&area);
    HBRUSH brush=CreateSolidBrush(dsl::theme::background);
    FillRect(dc,&area,brush); DeleteObject(brush); SetBkMode(dc,TRANSPARENT);
    if(solverActive_)solver_.paintEditorFrame(dc);
    SelectObject(dc,fonts_.title); SetTextColor(dc,dsl::theme::text);
    const wchar_t* title=L"DynSysLab";
    TextOutW(dc,scalePixels(window_,16),scalePixels(window_,10),title,lstrlenW(title));
    HPEN line=CreatePen(PS_SOLID,1,dsl::theme::border); auto old=SelectObject(dc,line);
    const int header=dsl::theme::metrics(window_).headerHeight;
    MoveToEx(dc,0,header-1,nullptr); LineTo(dc,area.right,header-1);
    SelectObject(dc,old); DeleteObject(line);
}

}  // namespace dsl::compute
